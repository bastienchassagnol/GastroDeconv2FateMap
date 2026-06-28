# Main functions

The Main function of CIBERSORT

## Usage

``` r
cibersort(sig_matrix, mixture_file, maxSize = 500, QN = TRUE)
```

## Arguments

- sig_matrix:

  sig_matrix file path to gene expression from isolated cells, or a
  matrix of expression profile of cells.

- mixture_file:

  mixture_file file path to heterogenous mixed expression file, or a
  matrix of heterogenous mixed expression

- maxSize:

  maximum size for the computation, to be passed to the
  future.global.maxSize. Default to 500 (in MB)

- QN:

  Perform quantile normalization or not (TRUE/FALSE)

## Examples

``` r
if (FALSE) { # \dontrun{
  ## example 1
  sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
  mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
  results <- cibersort(sig_matrix, mixture_file)
  ## example 2
  data(LM22)
  data(mixed_expr)
  results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr)
} # }
```
