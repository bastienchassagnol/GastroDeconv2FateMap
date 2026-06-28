# Pick a random anchor barcode for a phenotype

Samples one barcode identifier among organoids carrying the target
phenotype. Used to tag simulated samples and, for log-normal simulation,
to select a technical library-scaling factor.

## Usage

``` r
.pick_anchor_barcode(meta, phenotype_col, barcode_col, phenotype)
```

## Arguments

- meta:

  Cell metadata.

- phenotype_col:

  Phenotype column name.

- barcode_col:

  Barcode column name.

- phenotype:

  Target phenotype level.

## Value

Character scalar (barcode id).
