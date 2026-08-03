# Parallelism and memory

This document covers the CPU and memory settings of Prodes in full. **You do not need to read it to use Prodes.** The defaults are sensible, they are applied automatically, and a typical protein on a typical machine will run fine without you changing anything.

Read this if one of the following applies:

- You are running Prodes on a Linux server and want to control how many CPU cores it takes.
- You are processing hundreds or thousands of structures and want maximum throughput.
- You are calling Prodes from inside a web service, a job queue, or another program.
- You hit a memory error, or Prodes is using more RAM than you expected.

## The two kinds of parallelism

These are different things and it is worth keeping them apart, because combining them carelessly is the main way to overload a machine.

**One protein across many cores.** A single Prodes run splits the SASA, surface grid and shell phases across several worker processes. This is what makes one large protein finish quickly. It is controlled by `PRODES_N_WORKERS` and is the subject of most of this document.

**Many proteins at the same time, one core each.** You start several Prodes processes, each working on a different structure. This is what gives the best total throughput when you have a folder full of PDB files. Prodes does not do this for you; you do it with `multiprocessing`, GNU `parallel`, a shell loop, or a batch scheduler such as SLURM.

**Do not combine the two without doing the arithmetic.** Workers are cumulative. Ten processes at eight workers each will ask the machine for eighty cores. If your outer loop is already keeping the machine busy, set `PRODES_N_WORKERS=1` so each Prodes run stays on one core.

## Platform support

| Platform | Multi-core within one protein | Notes |
| --- | --- | --- |
| Linux | Yes, on by default | Uses half the logical cores unless told otherwise |
| Windows | No, always serial | The NumPy vectorisation speedup still applies |
| macOS | No, always serial | The NumPy vectorisation speedup still applies |

The workers share the parsed structure, the grid and the atom objects with their children through the `fork` start method. Windows and macOS do not provide it, so on those platforms Prodes runs serially. Nothing breaks and no setting is needed; the run is simply single-core.

This is not a limitation you can configure around, and it is why the benchmark headline of ~150x applies to Linux while Windows and macOS see ~30x.

## Choosing the number of workers

The default is `os.cpu_count() // 2`. `os.cpu_count()` reports *logical* cores, so on a hyperthreaded machine half of that is roughly the physical core count, which is the right target for numerical work.

Three ways to change it, in order of precedence:

1. **CLI flag**, which wins over everything else:

   ```bash
   python -m prodes in.pdb out.csv --n-workers 4
   ```

2. **Environment variable**, in your shell or in a `.env` file (see `.env.example`):

   ```bash
   export PRODES_N_WORKERS=4
   ```

3. **Nothing**, in which case you get half the logical cores.

`PRODES_N_WORKERS=1` forces a serial run even on Linux. Anything other than a positive integer raises a `ValueError` at startup rather than part-way through, so a typo fails in the first second.

There is also `--chunksize` / `PRODES_CHUNKSIZE`, the number of tasks dispatched to a worker at a time. The default of 1 gives the best load balancing, which matters because grid cells and atoms vary widely in how expensive they are. Raise it only if you have benchmarked your own workload and found that dispatch overhead dominates.

## What the parallelism actually buys you

On 1GPB (823 residues), the largest structure in the benchmark set, the calculation drops from **5146 s (86 minutes) to 30 s**, a **171x speedup**, in two independent steps:

- **28x from the NumPy vectorisation**, on a single core. This is the same work done without Python loops. It needs no extra hardware and applies on every platform.
- **A further ~6x from eight workers**, that is, eight CPU cores working on the one protein.

Eight cores give ~6x rather than 8x because two of the six pipeline phases are effectively serial, and the parallel phases pay the cost of starting a worker pool. 70% of the ideal is what this workload returns.

The full benchmark, including the data, the other structures and how the speedup varies with protein size, is in [docs/benchmark/benchmark_summary.md](benchmark/benchmark_summary.md).

## One process per protein, not one thread per protein

`calculate()` hands its work to the workers through module-level state. Two threads calling it in the same process would overwrite each other and return **wrong feature values rather than raising an error**. This is the one genuinely dangerous way to misuse Prodes.

Processes are unaffected, and every form of parallelism described in this document uses processes. But if you are calling Prodes from a web request handler, a `ThreadPoolExecutor`, or anything else that might call it from two threads at once, either serialise the calls or give each one its own process.

## Memory

The shell phase multiplies every charged atom against every surface point, and that intermediate array is by far the largest thing Prodes allocates. It is processed in chunks to stay within `PRODES_MEM_LIMIT_MB`, which is the budget for **one whole Prodes run** (default: 2048 MB) and is divided among the workers. Eight workers with the default get 256 MB each, so the run stays under 2 GB whatever the worker count.

**The budget is a ceiling, not an allocation.** A structure whose arrays are smaller than the budget runs in one chunk and uses only what it needs, which is why raising the setting on a small protein changes nothing at all. What a structure actually needs is roughly:

```
charged atoms x surface points x 48 bytes
```

Measured on this fork, single process, full feature set:

| Structure | Residues | Charged atoms | Surface points | Arrays needed | Peak process RAM |
| --- | --- | --- | --- | --- | --- |
| 1GDW | 130 | 68 | 7 039 | 22 MB | ~130 MB |
| 1GPB | 823 | 452 | 32 942 | 682 MB | ~920 MB |

Both counts grow with the protein, so the array grows with roughly the square of it: a complex twice the size of 1GPB needs something like 2.7 GB before chunking.

**Chunking is close to free.** The same 1GPB run took 127 s with a 2048 MB budget and 117 s with 64 MB, while peak RAM fell from 923 MB to 327 MB. There is therefore little to gain from a large budget and little to fear from a small one. Start low and raise it only if profiling says to.

Suggested settings, assuming the default of 8 workers per run:

| Machine | `PRODES_MEM_LIMIT_MB` | Peak RAM for the box |
| --- | --- | --- |
| 16 GB desktop, 1 protein at a time | 2048 (256 MB/worker) | ~3 GB |
| 128 GB server, 10 proteins at a time | 4096 (512 MB/worker) | ~50 GB |
| Under memory pressure, any machine | 512, or fewer workers | proportionally less |

The budget covers the chunked arrays only. The parsed structure, the grids and the surface points sit on top of it: about 170 MB on 1GPB, per Prodes process. The workers inherit those through `fork` and share the pages copy-on-write, so they cost rather less than a copy each, but not nothing. That overhead is why the server row above is ~50 GB rather than the 40 GB of budget alone.

### Setting the memory budget

- **Environment variable**, in your shell or in a `.env` file: `PRODES_MEM_LIMIT_MB=2048`
- **CLI flag**, which takes precedence: `python -m prodes in.pdb out.csv --mem-limit 2048`
- **Function argument**, when calling Prodes as a library: pass `mem_limit_mb` directly to `prodes.run.calculate(...)`. This is the recommended route inside worker processes, so that you are not relying on environment variables being inherited.

When running several proteins side by side, one process each, RAM is **cumulative** across those processes: each one gets its own budget. Divide the machine's memory by the number of concurrent runs, and leave headroom for the operating system.

## BLAS threads

Workers force single-threaded BLAS on themselves (`OPENBLAS_NUM_THREADS=1` and friends, applied via `threadpoolctl`), and the parent is held to one thread for as long as the workers exist. This is deliberate: parallelism should happen across processes rather than within them, and letting BLAS spawn its own threads inside each worker would oversubscribe the machine badly.

You only need to set these variables yourself if you also want to limit BLAS in the parent process outside the parallel phases.
