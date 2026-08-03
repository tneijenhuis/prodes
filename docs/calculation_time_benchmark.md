# Calculation Time Benchmark

## Overview

Surface property calculation time scales **exponentially** with the total
sequence length of all chains in the modelled structure. Users working with
large proteins or multimers should be aware that calculation times can become
substantial, and may wish to split large structures into smaller units
(individual chains or domains) where biologically meaningful.

![Surface property calculation time vs protein size](surface_calc_time_vs_length.png)

## Experiment

### Server specification

The benchmark was conducted on a Linux server with **125 GB of RAM** and a
single CPU core allocated to the calculation. No GPU acceleration was used.

### Dataset

Fifty-one *E. coli* proteins with homo-oligomeric structures were modelled
with [Boltz-2](https://github.com/jwohlwend/boltz) and used as input to
Prodes for surface property calculation. The total modelled sequence length
(sum of all chains in the multimer) ranged from **218 to 1,788 residues**.

### Software configuration

All calculations were performed with the NumPy-vectorised implementation
(see the [main README](../README.rst) for details on the vectorisation
work). Multiprocessing was **not** enabled — each calculation ran on a
single core.

### Results

The relationship between total sequence length and calculation time follows
an exponential curve:

> **y = 7.505 · exp(0.001988 · x) − 4.25**  &nbsp;&nbsp; (R² = 0.95)

where *x* is the total sequence length of all chains (amino acids) and *y*
is the calculation time in minutes.

### Interpretation

- **Small proteins (< 500 aa total):** calculations typically complete in
  under 15 minutes.
- **Medium proteins (500–1,000 aa total):** calculations range from ~15 to
  ~60 minutes.
- **Large proteins (> 1,000 aa total):** calculation time rises sharply,
  exceeding 2–4 hours for structures above ~1,500 residues.

The absolute times shown on the y-axis reflect single-core performance
without multiprocessing. Future commits enabling multiprocessing will reduce
these wall-clock times, but the **exponential scaling relationship** itself
is an inherent property of the algorithm and will persist regardless of
parallelisation strategy.

### Data

The raw timing data, including per-protein sequence lengths and oligomeric
states, is available in
[`surface_calc_times_with_length.csv`](surface_calc_times_with_length.csv).

A histogram of the calculation time distribution is also available:
[`surface_calc_time_histogram.png`](surface_calc_time_histogram.png).
