import logging
import multiprocessing as mp
import os
import sys
from collections.abc import Callable, Iterator
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from typing import Any

from threadpoolctl import threadpool_limits

logger = logging.getLogger(__name__)

# Environment variables read by the various BLAS backends when they start their
# thread pools. Set to "1" in every worker so that parallelism happens across
# processes rather than within them.
BLAS_THREAD_ENV_VARS = (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def parallelism_enabled() -> bool:
    """Return True when parallel processing is available on this platform.

    The worker design relies on the fork start method to share the structure,
    grid and atom objects with the children, so only Linux is supported.
    Windows and macOS always use the serial path.
    """

    return sys.platform.startswith("linux")


def requested_n_workers() -> int:
    """Return the worker count asked for via PRODES_N_WORKERS, ignoring the platform.

    Defaults to half the logical CPU count, which approximates the physical core
    count on hyperthreaded systems (os.cpu_count reports logical cores). Kept
    separate from n_workers so that a malformed setting is reported on every
    platform rather than only where parallelism is available.

    Raises:
        ValueError: if PRODES_N_WORKERS is set to anything but a positive integer.
    """

    setting = os.getenv("PRODES_N_WORKERS")
    if setting is None or not setting.strip():
        return max(1, (os.cpu_count() or 2) // 2)

    try:
        workers = int(setting)
    except ValueError:
        raise ValueError(f"PRODES_N_WORKERS must be a positive integer, got {setting!r}") from None

    if workers < 1:
        raise ValueError(f"PRODES_N_WORKERS must be >= 1, got {workers}")

    return workers


def n_workers() -> int:
    """Return the number of worker processes to use for this run.

    Returns 1 on platforms without parallelism support and when
    PRODES_N_WORKERS=1, either of which forces the serial path.
    """

    workers = requested_n_workers()

    return workers if parallelism_enabled() else 1


def chunksize() -> int:
    """Return the chunk size for ProcessPoolExecutor.map calls.

    Defaults to 1, which gives the best load balancing when the cost of
    individual tasks varies widely, as it does for grid cells and atoms.
    Raise it via PRODES_CHUNKSIZE only if benchmarking shows that dispatch
    overhead dominates.

    Raises:
        ValueError: if PRODES_CHUNKSIZE is set to anything but a positive integer.
    """

    setting = os.getenv("PRODES_CHUNKSIZE")
    if setting is None or not setting.strip():
        return 1

    try:
        size = int(setting)
    except ValueError:
        raise ValueError(f"PRODES_CHUNKSIZE must be a positive integer, got {setting!r}") from None

    if size < 1:
        raise ValueError(f"PRODES_CHUNKSIZE must be >= 1, got {size}")

    return size


def run_mem_limit_mb(setting: float | None = None) -> float:
    """Return the transient-array budget for one whole Prodes run, in MB.

    This is the budget for the run as a whole, not for one worker: see
    worker_mem_limit_mb, which is what the calculations actually get.

    Args:
        setting: an explicit budget, which wins over the environment variable.
            Pass None to read PRODES_MEM_LIMIT_MB (2048 MB if unset).

    Raises:
        ValueError: if the budget is not a positive number.
    """

    if setting is None:
        raw = os.getenv("PRODES_MEM_LIMIT_MB")
        if raw is None or not raw.strip():
            return 2048.0
        try:
            setting = float(raw)
        except ValueError:
            raise ValueError(f"PRODES_MEM_LIMIT_MB must be a positive number, got {raw!r}") from None

    setting = float(setting)
    if setting <= 0:
        raise ValueError(f"PRODES_MEM_LIMIT_MB must be > 0, got {setting}")

    return setting


def worker_mem_limit_mb(setting: float | None = None) -> float:
    """Return the share of the run's memory budget that one worker may use.

    Every worker runs the same chunked calculations, so a budget applied per
    worker would be multiplied by the worker count in real memory. Dividing it
    here keeps PRODES_MEM_LIMIT_MB meaning what it says: the ceiling on transient
    arrays for the run, whatever the worker count. A smaller share means smaller
    chunks, which costs time rather than correctness.

    Args:
        setting: an explicit whole-run budget in MB, or None to read the
            environment variable.
    """

    return max(1.0, run_mem_limit_mb(setting) / n_workers())


def validate_settings(mem_limit_mb: float | None = None) -> None:
    """Check every parallelism setting once, so a bad value fails immediately.

    The settings are otherwise read at the point of use, which is partway
    through the first phase, after parsing and grid construction have been paid
    for. Call this at startup.

    Args:
        mem_limit_mb: an explicit memory budget to validate instead of the
            environment variable, as passed to calculate() or --mem-limit.

    Raises:
        ValueError: if PRODES_N_WORKERS, PRODES_CHUNKSIZE or the memory budget is
            set to anything other than a valid value.
    """

    requested_n_workers()
    chunksize()
    run_mem_limit_mb(mem_limit_mb)


def limit_threads() -> None:
    """Force single-threaded BLAS in the calling process.

    Used as the worker pool initializer. Forking a parent whose BLAS backend
    already has a live thread pool leaves that pool broken in the child, which
    can deadlock. The environment variables cover backends that read them when
    they start up, and threadpool_limits additionally reconfigures pools that
    are already running, which is the case immediately after a fork.
    """

    for variable in BLAS_THREAD_ENV_VARS:
        os.environ[variable] = "1"

    threadpool_limits(limits=1)


@contextmanager
def worker_pool(max_workers: int) -> Iterator[ProcessPoolExecutor | None]:
    """Yield a worker pool, or None when the caller should run serially.

    None is yielded when parallelism is unavailable on this platform, when a
    single worker was requested, and when constructing the pool itself raises.
    The last case is logged rather than raised, so that an unusable pool degrades
    to serial execution instead of failing the run.

    Constructing a ProcessPoolExecutor does not start any process, so a fork that
    fails does not surface here: it surfaces at the first map call, where
    run_tasks catches it and falls back to serial. Both routes end up serial, but
    they are two different fallbacks and only the first one is this function's.

    Args:
        max_workers: number of worker processes to start.
    """

    if max_workers <= 1 or not parallelism_enabled():
        yield None
        return

    try:
        pool = ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=mp.get_context("fork"),
            initializer=limit_threads,
        )
    except Exception as error:
        logger.warning("Could not create a worker pool, falling back to serial execution: %s", error)
        yield None
        return

    try:
        yield pool
    finally:
        pool.shutdown(wait=True)


def run_tasks(worker: Callable[[int], Any], n_tasks: int, description: str) -> list:
    """Run worker(index) for every index in range(n_tasks), in parallel where possible.

    The worker must be a module-level function taking a single integer index and
    reading everything else from module-level state that the parent populated
    before the pool was created. Fork gives the children that state through
    copy-on-write, so nothing but the index and the result crosses the pipe.
    Any failure in the parallel path is logged and the whole run is repeated
    serially, since a partially consumed map leaves no usable results.

    The parent's own BLAS pool is capped to one thread for as long as the workers
    exist. The deadlock this design has to avoid comes from forking a process
    whose BLAS backend has live threads, and by the time the first pool is
    created the parent has done a great deal of NumPy work. Capping the child in
    the pool initializer repairs the pool after the fork; capping the parent here
    closes the window instead. The cap is a context manager, so the parent's
    original thread count is restored before the serial path runs.

    Args:
        worker: module-level function mapping a task index to a result.
        n_tasks: number of tasks to run; the indices are range(n_tasks).
        description: short phrase naming the work, used in the fallback warning.
    """

    indices = range(n_tasks)
    workers = n_workers()

    # n_workers is already 1 wherever parallelism is unavailable, so this also
    # keeps the parent's thread count untouched on the serial platforms.
    if workers > 1:
        with threadpool_limits(limits=1):
            with worker_pool(workers) as pool:
                if pool is not None:
                    try:
                        return list(pool.map(worker, indices, chunksize=chunksize()))
                    except Exception as error:
                        logger.warning("Parallel %s failed, falling back to serial execution: %s", description, error)

    return [worker(index) for index in indices]
