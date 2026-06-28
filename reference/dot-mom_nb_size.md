# Per-gene NB dispersion via method of moments

Computes the negative-binomial `size` (dispersion) parameter from a
vector of counts using method-of-moments, as in HADACA3
pseudo-single-cell generation.

## Usage

``` r
.mom_nb_size(x)
```

## Arguments

- x:

  Numeric vector of counts for one gene.

## Value

Dispersion `size` (theta); may be `Inf`.

## Details

Method-of-moments negative-binomial size from counts

For counts with mean \\m\\ and variance \\v\\, under \\X \sim
\mathrm{NB}(\mu=m, \mathrm{size}=\theta)\\: \$\$\mathrm{Var}(X) = m +
\frac{m^2}{\theta} \quad\Rightarrow\quad \hat{\theta} = \frac{m^2}{v -
m}\$\$ when \\v \> m\\. Otherwise \\\hat{\theta} = \infty\\ (Poisson
limit).

## References

HADACA3 in silico pseudo-bulk workflow
(<https://github.com/bioinfo-LIG/hadaca3_framework>), which fits
[`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html) per gene;
method-of-moments is used here for speed at scale.
