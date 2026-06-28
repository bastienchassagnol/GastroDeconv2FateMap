# Normalise a composition vector to the unit simplex

Rescales a non-negative numeric vector so that its entries sum to 1
(`\\sum_j x_j = 1`), i.e. the unit simplex constraint for compositions.

## Usage

``` r
normalise_cell_estimates(x)
```

## Arguments

- x:

  Numeric vector containing non-negative cell estimates.

## Value

A numeric vector on the unit simplex with the same names as `x`.

## Details

If the input already satisfies the simplex constraint (within
`sqrt(.Machine$double.eps)`), it is returned unchanged. Otherwise, it is
divided by its sum.
