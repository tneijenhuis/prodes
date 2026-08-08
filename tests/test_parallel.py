"""Tests for the parallel execution infrastructure and its serial fallbacks.

The parallel path itself is Linux only, because the workers rely on the fork
start method to share the structure, grid and atom objects. Everything that does
not need a real fork is tested on every platform by monkeypatching
``parallelism_enabled``; the tests that genuinely need worker processes are
skipped elsewhere and run in CI, which is Linux.
"""

import concurrent.futures
import logging
import os
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pytest
from threadpoolctl import threadpool_info, threadpool_limits

from prodes import parallel
from prodes.calculations import sasa
from prodes.run import AVERAGE_CHARGE_GRID_STATE, SHELL_STATE, SURFACE_GRID_STATE, calculate

on_linux_only = pytest.mark.skipif(not sys.platform.startswith("linux"), reason="the parallel path requires the fork start method, which is Linux only")

PDB_PATH = "tests/data/ARH96693.pdb.zip"


# -- Worker count --


def test_requested_n_workers_defaults_to_half_the_logical_cpus(monkeypatch):
    """Without PRODES_N_WORKERS the default is half the logical CPU count."""
    monkeypatch.delenv("PRODES_N_WORKERS", raising=False)
    assert parallel.requested_n_workers() == max(1, (os.cpu_count() or 2) // 2)


def test_requested_n_workers_reads_the_env_var(monkeypatch):
    """An explicit PRODES_N_WORKERS overrides the default."""
    monkeypatch.setenv("PRODES_N_WORKERS", "4")
    assert parallel.requested_n_workers() == 4


def test_requested_n_workers_ignores_a_blank_setting(monkeypatch):
    """An empty PRODES_N_WORKERS falls back to the default rather than raising."""
    monkeypatch.setenv("PRODES_N_WORKERS", "   ")
    assert parallel.requested_n_workers() == max(1, (os.cpu_count() or 2) // 2)


def test_requested_n_workers_rejects_a_non_integer(monkeypatch):
    """A non-integer PRODES_N_WORKERS raises a ValueError naming the variable and the value."""
    monkeypatch.setenv("PRODES_N_WORKERS", "abc")
    with pytest.raises(ValueError, match="PRODES_N_WORKERS must be a positive integer, got 'abc'"):
        parallel.requested_n_workers()


@pytest.mark.parametrize("setting", ["0", "-1"])
def test_requested_n_workers_rejects_a_non_positive_integer(monkeypatch, setting):
    """Zero or a negative PRODES_N_WORKERS raises a ValueError."""
    monkeypatch.setenv("PRODES_N_WORKERS", setting)
    with pytest.raises(ValueError, match="PRODES_N_WORKERS must be >= 1"):
        parallel.requested_n_workers()


def test_n_workers_is_one_when_parallelism_is_unavailable(monkeypatch):
    """On a platform without parallelism support the worker count collapses to 1."""
    monkeypatch.setenv("PRODES_N_WORKERS", "8")
    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: False)
    assert parallel.n_workers() == 1


def test_n_workers_honours_the_env_var_where_parallelism_is_available(monkeypatch):
    """Where parallelism is available the requested worker count is used as is."""
    monkeypatch.setenv("PRODES_N_WORKERS", "8")
    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: True)
    assert parallel.n_workers() == 8


def test_parallelism_is_enabled_only_on_linux():
    """Parallelism is gated on Linux, since the workers need the fork start method."""
    assert parallel.parallelism_enabled() == sys.platform.startswith("linux")


# -- Chunk size --


def test_chunksize_defaults_to_one(monkeypatch):
    """Without PRODES_CHUNKSIZE tasks are dispatched one at a time for load balancing."""
    monkeypatch.delenv("PRODES_CHUNKSIZE", raising=False)
    assert parallel.chunksize() == 1


def test_chunksize_reads_the_env_var(monkeypatch):
    """An explicit PRODES_CHUNKSIZE overrides the default."""
    monkeypatch.setenv("PRODES_CHUNKSIZE", "5")
    assert parallel.chunksize() == 5


def test_chunksize_rejects_a_non_integer(monkeypatch):
    """A non-integer PRODES_CHUNKSIZE raises a ValueError naming the variable and the value."""
    monkeypatch.setenv("PRODES_CHUNKSIZE", "big")
    with pytest.raises(ValueError, match="PRODES_CHUNKSIZE must be a positive integer, got 'big'"):
        parallel.chunksize()


def test_chunksize_rejects_a_non_positive_integer(monkeypatch):
    """A PRODES_CHUNKSIZE below 1 raises a ValueError."""
    monkeypatch.setenv("PRODES_CHUNKSIZE", "0")
    with pytest.raises(ValueError, match="PRODES_CHUNKSIZE must be >= 1"):
        parallel.chunksize()


# -- Memory budget --


def test_run_mem_limit_defaults_to_2048(monkeypatch):
    """Without PRODES_MEM_LIMIT_MB the whole-run budget is 2048 MB."""
    monkeypatch.delenv("PRODES_MEM_LIMIT_MB", raising=False)
    assert parallel.run_mem_limit_mb() == 2048.0


def test_run_mem_limit_reads_the_env_var(monkeypatch):
    """An explicit PRODES_MEM_LIMIT_MB overrides the default."""
    monkeypatch.setenv("PRODES_MEM_LIMIT_MB", "512")
    assert parallel.run_mem_limit_mb() == 512.0


def test_run_mem_limit_prefers_an_explicit_setting(monkeypatch):
    """A budget passed in code wins over the environment, as --mem-limit does."""
    monkeypatch.setenv("PRODES_MEM_LIMIT_MB", "512")
    assert parallel.run_mem_limit_mb(256) == 256.0


@pytest.mark.parametrize("setting", ["abc", "0", "-1"])
def test_run_mem_limit_rejects_a_bad_setting(monkeypatch, setting):
    """A non-numeric or non-positive budget raises a ValueError naming the variable."""
    monkeypatch.setenv("PRODES_MEM_LIMIT_MB", setting)
    with pytest.raises(ValueError, match="PRODES_MEM_LIMIT_MB must be"):
        parallel.run_mem_limit_mb()


def test_worker_mem_limit_divides_the_budget_by_the_worker_count(monkeypatch):
    """Every worker runs the same chunked calculation, so each gets a share of the budget.

    Handing the whole budget to each worker is what makes peak memory
    n_workers times what the setting says.
    """
    monkeypatch.setattr(parallel, "n_workers", lambda: 8)
    assert parallel.worker_mem_limit_mb(4096) == 512.0


def test_worker_mem_limit_is_the_whole_budget_when_serial(monkeypatch):
    """A serial run has nothing to share the budget with."""
    monkeypatch.setattr(parallel, "n_workers", lambda: 1)
    assert parallel.worker_mem_limit_mb(4096) == 4096.0


def test_worker_mem_limit_never_reaches_zero(monkeypatch):
    """A tiny budget spread over many workers still leaves each of them something."""
    monkeypatch.setattr(parallel, "n_workers", lambda: 64)
    assert parallel.worker_mem_limit_mb(8) == 1.0


# -- Startup validation --


def test_validate_settings_accepts_a_sane_environment(monkeypatch):
    """Valid settings pass quietly."""
    monkeypatch.setenv("PRODES_N_WORKERS", "4")
    monkeypatch.setenv("PRODES_CHUNKSIZE", "2")
    monkeypatch.setenv("PRODES_MEM_LIMIT_MB", "1024")

    parallel.validate_settings()


@pytest.mark.parametrize(
    ("variable", "value"),
    [("PRODES_N_WORKERS", "abc"), ("PRODES_CHUNKSIZE", "0"), ("PRODES_MEM_LIMIT_MB", "-5")],
)
def test_validate_settings_rejects_a_bad_value(monkeypatch, variable, value):
    """Each setting is checked up front rather than partway through the first phase."""
    monkeypatch.setenv(variable, value)
    with pytest.raises(ValueError, match=variable):
        parallel.validate_settings()


def test_calculate_validates_before_doing_any_work(monkeypatch, tmp_path):
    """A bad setting fails immediately, not after the structure has been parsed."""

    monkeypatch.setenv("PRODES_N_WORKERS", "abc")

    def fail_if_called(*args, **kwargs):
        raise AssertionError("the pipeline started before the settings were checked")

    monkeypatch.setattr("prodes.run.prepare_structure", fail_if_called)

    with pytest.raises(ValueError, match="PRODES_N_WORKERS"):
        calculate(PDB_PATH, str(tmp_path / "output.csv"))


# -- BLAS thread limiting --


@pytest.fixture
def restore_blas_threads():
    """Puts the process's BLAS thread counts back after a test calls limit_threads.

    limit_threads is meant to run in a worker that is about to do its own work
    and then exit. Called in the test process it has no teardown of its own, so
    every later test in the session would run against single-threaded BLAS.
    """

    original = threadpool_limits()
    try:
        yield
    finally:
        original.restore_original_limits()


def test_limit_threads_sets_every_blas_variable(monkeypatch, restore_blas_threads):
    """limit_threads pins all known BLAS thread variables to 1."""
    for variable in parallel.BLAS_THREAD_ENV_VARS:
        monkeypatch.delenv(variable, raising=False)

    parallel.limit_threads()

    assert all(os.environ[variable] == "1" for variable in parallel.BLAS_THREAD_ENV_VARS)


def test_limit_threads_overrides_an_existing_setting(monkeypatch, restore_blas_threads):
    """A thread count inherited from the parent is overridden rather than kept."""
    monkeypatch.setenv("OMP_NUM_THREADS", "8")

    parallel.limit_threads()

    assert os.environ["OMP_NUM_THREADS"] == "1"


# -- Worker pool --


def test_worker_pool_yields_none_for_a_single_worker(monkeypatch):
    """Asking for one worker means serial execution, signalled by yielding None."""
    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: True)
    with parallel.worker_pool(1) as pool:
        assert pool is None


def test_worker_pool_yields_none_when_parallelism_is_unavailable(monkeypatch):
    """On a platform without parallelism support no pool is created."""
    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: False)
    with parallel.worker_pool(8) as pool:
        assert pool is None


def test_worker_pool_degrades_to_serial_when_creation_fails(monkeypatch, caplog):
    """A pool that cannot be created is logged and degrades to serial rather than raising."""

    def refuse_to_start(*args, **kwargs):
        raise OSError("no worker processes available")

    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: True)
    # Off Linux the fork context is itself unavailable, which would fail before
    # the executor is ever reached, so it is stubbed out to isolate the fallback.
    monkeypatch.setattr(parallel.mp, "get_context", lambda method: None)
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", refuse_to_start)

    with caplog.at_level(logging.WARNING, logger="prodes.parallel"):
        with parallel.worker_pool(4) as pool:
            assert pool is None

    assert "falling back to serial" in caplog.text
    assert "no worker processes available" in caplog.text


def test_worker_pool_degrades_to_serial_without_a_fork_context(monkeypatch, caplog):
    """A platform that claims parallelism but cannot fork still degrades rather than raising."""

    def no_such_context(method):
        raise ValueError(f"cannot find context for {method!r}")

    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: True)
    monkeypatch.setattr(parallel.mp, "get_context", no_such_context)

    with caplog.at_level(logging.WARNING, logger="prodes.parallel"):
        with parallel.worker_pool(4) as pool:
            assert pool is None

    assert "falling back to serial" in caplog.text


# -- Task dispatch --


def square(index: int) -> int:
    """Module-level worker used to exercise run_tasks without a closure."""
    return index * index


def test_run_tasks_runs_every_index(monkeypatch):
    """run_tasks covers range(n_tasks) and preserves the input order."""
    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: False)
    assert parallel.run_tasks(square, 5, "squares") == [0, 1, 4, 9, 16]


def test_run_tasks_returns_an_empty_list_for_no_tasks(monkeypatch):
    """A calculation with nothing to do returns an empty list rather than failing."""
    monkeypatch.setattr(parallel, "parallelism_enabled", lambda: False)
    assert parallel.run_tasks(square, 0, "squares") == []


def blas_thread_counts():
    """Thread counts of every BLAS pool threadpoolctl can see in this process."""

    return [pool["num_threads"] for pool in threadpool_info()]


def test_run_tasks_pins_the_parents_blas_threads_while_workers_exist(monkeypatch):
    """The parent forks single-threaded and gets its own thread count back afterwards.

    Forking a parent whose BLAS backend has live threads is what deadlocks, and
    the pool initializer only repairs the child after the fork has happened.
    """

    before = blas_thread_counts()
    if not before:
        pytest.skip("no BLAS backend visible to threadpoolctl in this environment")

    seen = {}

    @contextmanager
    def recording_worker_pool(max_workers):
        seen["during"] = blas_thread_counts()
        yield None

    monkeypatch.setattr(parallel, "n_workers", lambda: 4)
    monkeypatch.setattr(parallel, "worker_pool", recording_worker_pool)

    assert parallel.run_tasks(square, 3, "squares") == [0, 1, 4]

    assert seen["during"] == [1] * len(before)
    assert blas_thread_counts() == before


def test_run_tasks_leaves_the_parents_blas_threads_alone_when_serial(monkeypatch):
    """A serial run never forks, so there is nothing to protect and nothing to change."""

    before = blas_thread_counts()
    if not before:
        pytest.skip("no BLAS backend visible to threadpoolctl in this environment")

    monkeypatch.setattr(parallel, "n_workers", lambda: 1)
    assert parallel.run_tasks(square, 3, "squares") == [0, 1, 4]

    assert blas_thread_counts() == before


def test_run_tasks_falls_back_to_serial_when_the_pool_fails(monkeypatch, caplog):
    """A pool that fails mid-map is logged and the whole run is repeated serially."""

    class BrokenPool:
        def map(self, *args, **kwargs):
            raise concurrent.futures.BrokenExecutor("worker died")

    @contextmanager
    def broken_worker_pool(max_workers):
        yield BrokenPool()

    monkeypatch.setattr(parallel, "n_workers", lambda: 4)
    monkeypatch.setattr(parallel, "worker_pool", broken_worker_pool)

    with caplog.at_level(logging.WARNING, logger="prodes.parallel"):
        results = parallel.run_tasks(square, 4, "squares")

    assert results == [0, 1, 4, 9]
    assert "Parallel squares failed, falling back to serial" in caplog.text


# -- Shared state --


def test_calculate_leaves_no_state_behind(tmp_path):
    """The state a phase hands to its workers is dropped as soon as the phase is done.

    The dicts hold the structure, both grids and the SASA task arrays, which on a
    large structure is the bulk of the run's memory. Keeping them would hold one
    structure alive until the next call replaced it, and the last one for the
    life of the process.
    """

    calculate(PDB_PATH, str(tmp_path / "output.csv"), full_features=True)

    assert SURFACE_GRID_STATE == {}
    assert AVERAGE_CHARGE_GRID_STATE == {}
    assert SHELL_STATE == {}
    assert sasa.SASA_TASKS == {}


# -- End to end equivalence, which needs real worker processes --


def features_for(monkeypatch, **env):
    """Runs the full pipeline under the given environment and returns the feature row.

    Args:
        monkeypatch: pytest fixture used to scope the environment to one test.
        **env: environment variables to set for this run.
    """

    import pandas as pd

    for name, value in env.items():
        monkeypatch.setenv(name, value)

    with tempfile.TemporaryDirectory() as directory:
        out_file = str(Path(directory) / "output.csv")
        calculate(PDB_PATH, out_file, full_features=True)
        return pd.read_csv(out_file).iloc[0].drop("ID").astype(float)


def assert_features_match(serial, parallel_run):
    """Asserts that two feature rows agree to within floating point tolerance."""

    assert set(serial.index) == set(parallel_run.index)
    for column in serial.index:
        assert parallel_run[column] == pytest.approx(serial[column], rel=1e-6, abs=1e-9), f"{column}: serial={serial[column]}, parallel={parallel_run[column]}"


@on_linux_only
def test_parallel_matches_serial(monkeypatch):
    """The parallel path reproduces the serial feature values."""
    serial = features_for(monkeypatch, PRODES_N_WORKERS="1")
    parallel_run = features_for(monkeypatch, PRODES_N_WORKERS="4")
    assert_features_match(serial, parallel_run)


@on_linux_only
def test_chunksize_does_not_change_the_result(monkeypatch):
    """Dispatching several tasks per chunk gives the same values as one at a time."""
    single = features_for(monkeypatch, PRODES_N_WORKERS="4", PRODES_CHUNKSIZE="1")
    batched = features_for(monkeypatch, PRODES_N_WORKERS="4", PRODES_CHUNKSIZE="5")
    assert_features_match(single, batched)


@on_linux_only
def test_worker_pool_starts_real_processes(monkeypatch):
    """The pool runs work in separate processes with BLAS pinned to one thread."""
    monkeypatch.delenv("PRODES_N_WORKERS", raising=False)

    with parallel.worker_pool(2) as pool:
        assert pool is not None
        results = list(pool.map(worker_fingerprint, range(4)))

    parent = os.getpid()
    assert all(pid != parent for pid, _ in results)
    assert {threads for _, threads in results} == {"1"}


def worker_fingerprint(index: int):
    """Module-level worker reporting its process id and BLAS thread setting."""
    return os.getpid(), os.environ.get("OMP_NUM_THREADS")


@on_linux_only
def test_repeated_pool_creation_does_not_deadlock(monkeypatch):
    """Creating and tearing down pools repeatedly must not hang after NumPy has used BLAS.

    Forking a parent whose BLAS backend already has live threads is the classic
    source of intermittent deadlocks, and it only shows up after several pools
    have been created. Kept deliberately cheap so it is affordable in CI.
    """

    monkeypatch.delenv("PRODES_N_WORKERS", raising=False)

    # Make sure the parent's BLAS thread pool is actually running before forking
    np.linalg.eigh(np.random.default_rng(0).random((200, 200)))

    for _ in range(10):
        assert parallel.run_tasks(square, 8, "squares") == [index * index for index in range(8)]
