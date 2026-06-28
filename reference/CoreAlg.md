# Core algorithm

The core algorithm of CIBERSORT which uses nu-SVR
([`e1071::svm`](https://rdrr.io/pkg/e1071/man/svm.html)).

## Usage

``` r
CoreAlg(X, y, maxSize = 500)
```

## Arguments

- X:

  Cell-specific gene expression matrix.

- y:

  Mixed expression vector for one sample.

- maxSize:

  Maximum global object size (MB) for parallel workers.

## Value

A list with `w` (weights), `mix_rmse`, and `mix_r`.
