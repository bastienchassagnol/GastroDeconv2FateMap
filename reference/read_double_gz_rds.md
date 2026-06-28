# Read a doubly gzip-compressed RDS file

Loads an R object from an `.rds.gz` file that has been compressed twice
with gzip. A single [`readRDS`](https://rdrr.io/r/base/readRDS.html)
call on [`gzfile`](https://rdrr.io/r/base/connections.html) only removes
the outer compression layer; the remaining bytes are still gzip-encoded,
which triggers an *unknown input format* error. This helper decompresses
both layers in memory, then calls
[`readRDS`](https://rdrr.io/r/base/readRDS.html).

Some GEO supplementary archives (e.g. per-sample files inside
`GSE250136_RAW.tar`) are stored in this double-gzip form.

## Usage

``` r
read_double_gz_rds(path)
```

## Arguments

- path:

  Character scalar. Path to the `.rds.gz` file.

## Value

The R object stored in the RDS file (class depends on the file; often a
Seurat object for single-cell supplementary data).

## Details

Peak memory use is roughly the size of the first decompressed gzip layer
plus the second. For large objects, ensure sufficient RAM is available.

## See also

[readRDS](https://rdrr.io/r/base/readRDS.html),
[gzfile](https://rdrr.io/r/base/connections.html),
[gzcon](https://rdrr.io/r/base/gzcon.html),
[rawConnection](https://rdrr.io/r/base/rawConnection.html),
[readBin](https://rdrr.io/r/base/readBin.html)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- read_double_gz_rds(
  "./data/raw/GSE250136/GSM7974412_df48_final.rds.gz"
)
} # }
```
