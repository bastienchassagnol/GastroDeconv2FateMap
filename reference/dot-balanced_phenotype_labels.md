# Balanced phenotype labels for simulated samples

Allocates `n_samples` across phenotype levels as evenly as possible.
When `n_samples` is not divisible by the number of phenotypes, the first
`n_samples %% K` levels receive one extra sample each.

## Usage

``` r
.balanced_phenotype_labels(n_samples, phenotypes)
```

## Arguments

- n_samples:

  Total number of pseudo-bulk samples to generate.

- phenotypes:

  Character vector of phenotype levels.

## Value

Character vector of length `n_samples` with (near-)balanced phenotype
labels.
