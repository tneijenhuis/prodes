# Relationship between protein size and calculation time

## Overview

For larger complexes, surface property calculation time scales **exponentially** with the total
sequence length of all chains in the modelled structure. Users working with
large proteins or multimers should be aware that calculation times can become
substantial, and may wish to split large structures into smaller units
(individual chains or domains) where biologically meaningful.

![Surface property calculation time vs protein size](surface_calc_time_vs_length.png)

## Experiment

### Server specification

The benchmark was conducted on a Linux server (Server1) with **125 GB of RAM** and a
single CPU core allocated to the calculation. No GPU acceleration was used.

### Dataset

Fifty-one *E. coli* proteins with homo-oligomeric structures were modelled
with [Boltz-2](https://github.com/jwohlwend/boltz) and used as input to
Prodes for surface property calculation. The total modelled sequence length
(sum of all chains in the multimer) ranged from **218 to 1,788 residues**.
Residue pKa values were assigned with propKa.

### Software configuration

All calculations were performed with the NumPy-vectorised implementation
(see the [main README](../README.rst) for details on the vectorisation
work). Multiprocessing was **not** enabled — each calculation ran on a
single core. Not all improvements from the vectorisation had been fully
implemented yet.

### Results

The relationship between total sequence length and calculation time follows
an exponential curve:

> **y = 7.505 · exp(0.001988 · x) − 4.25**  &nbsp;&nbsp; (R² = 0.95)

where *x* is the total sequence length of all chains (amino acids) and *y*
is the calculation time in minutes.

### Interpretation

The absolute times shown on the y-axis reflect single-core performance
without multiprocessing. Multiprocessing utilizing multiple CPU cores has since been added, and the wall-clock times below no longer reflect what the current code does: see
[`benchmark/benchmark_summary.md`](benchmark/benchmark_summary.md) for a direct comparison against
the original prodes, where the largest structure measured drops from 80+ minutes to ~0.5 minutes.

### Data

The raw timing data, including per-protein sequence lengths and oligomeric
states, is available in
[`surface_calc_times_with_length.csv`](surface_calc_times_with_length.csv).

A histogram of the calculation time distribution is also available:
[`surface_calc_time_histogram.png`](surface_calc_time_histogram.png).
