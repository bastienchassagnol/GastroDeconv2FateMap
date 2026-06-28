# Configure a furrr parallel backend for CIBERSORT

Configure a furrr parallel backend for CIBERSORT

## Usage

``` r
.enable_parallel_cibersort(nThreads = NULL, maxSize = 500)
```

## Arguments

- nThreads:

  Number of parallel workers. Defaults to `availableCores() - 2`.

- maxSize:

  Maximum global object size (MB) for future exports.
