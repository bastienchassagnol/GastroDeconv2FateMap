"""scCODA compositional analysis — Luque 120h TLS vs neural_bias (pertpy)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pertpy as pt
from anndata import read_h5ad

# %%
STUDY = "luque_GSE250136"
TIMEPOINT = "120h"
CONTRAST = "neural_bias vs TLS (reference morphotype: TLS)"
REFERENCE_CELL_TYPE = "pluripotent"

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "biological-exploration" / "pertpy"
OUT.mkdir(parents=True, exist_ok=True)

H5AD = ROOT / "data" / "intermediate" / "luque_120h_sccoda.h5ad"
FORMULA = "C(Morphotype, Treatment('TLS'))"
EST_FDR = 0.05
CI_PROB = 0.95

# %%
adata = read_h5ad(H5AD)
sccoda = pt.tl.Sccoda()

mdata = sccoda.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="luque_cluster_annotation",
    sample_identifier="Sample.barcode",
    covariate_obs=["Morphotype"],
)

# %%
sccoda.plot_stacked_barplot(
    mdata,
    feature_name="Morphotype",
    level_order=["TLS", "neural_bias"],
    figsize=(7, 4),
    return_fig=True,
)
plt.tight_layout()
plt.savefig(OUT / f"{STUDY}_{TIMEPOINT}_sccoda_stacked_barplot.pdf")
plt.close()

sccoda.plot_rel_abundance_dispersion_plot(
    mdata,
    abundant_threshold=0.95,
    label_cell_types=True,
    figsize=(7, 5),
    return_fig=True,
)
plt.tight_layout()
plt.savefig(OUT / f"{STUDY}_{TIMEPOINT}_sccoda_reference_dispersion.pdf")
plt.close()

# %%
mdata = sccoda.prepare(
    mdata,
    formula=FORMULA,
    reference_cell_type=REFERENCE_CELL_TYPE,
)

sccoda.run_nuts(
    mdata,
    num_samples=20_000,
    num_warmup=5_000,
    rng_key=42,
)

# ci_prob=0.95: pertpy/az.summary default is ~89% ETI if omitted
sccoda.set_fdr(mdata, est_fdr=EST_FDR, ci_prob=CI_PROB)

intercepts = sccoda.get_intercept_df(mdata).reset_index()
effects = sccoda.get_effect_df(mdata).reset_index()
credible = sccoda.credible_effects(mdata, est_fdr=EST_FDR).rename("credible")
credible = credible.reset_index()

eti_low = next(c for c in effects.columns if c.startswith("ETI") and "lower" in c)
eti_high = next(c for c in effects.columns if c.startswith("ETI") and "upper" in c)
log2 = np.log(2)
effects["relative_log2_ratio"] = effects["Final Parameter"] / log2
effects["relative_log2_low"] = effects[eti_low] / log2
effects["relative_log2_high"] = effects[eti_high] / log2

for tbl in (intercepts, effects, credible):
    tbl.insert(0, "study", STUDY)
    tbl.insert(1, "timepoint", TIMEPOINT)
    tbl.insert(2, "contrast", CONTRAST)
    tbl.insert(3, "reference_cell_type", REFERENCE_CELL_TYPE)
effects.insert(4, "ci_prob", CI_PROB)

intercepts.to_csv(OUT / f"{STUDY}_{TIMEPOINT}_sccoda_intercepts.csv", index=False)
effects.to_csv(OUT / f"{STUDY}_{TIMEPOINT}_sccoda_effects_extended.csv", index=False)
credible.to_csv(OUT / f"{STUDY}_{TIMEPOINT}_sccoda_credible_effects.csv", index=False)

print(sccoda.summary(mdata, extended=True, ci_prob=CI_PROB))

# %%
idata = sccoda.make_arviz(mdata)

diagnostics = az.summary(
    idata,
    var_names=["alpha", "beta"],
    ci_prob=CI_PROB,
)
diagnostics = diagnostics.reset_index(names="parameter")
diagnostics.insert(0, "study", STUDY)
diagnostics.insert(1, "timepoint", TIMEPOINT)
diagnostics.insert(2, "contrast", CONTRAST)
diagnostics.insert(3, "reference_cell_type", REFERENCE_CELL_TYPE)
diagnostics.insert(4, "ci_prob", CI_PROB)
diagnostics.to_csv(OUT / f"{STUDY}_{TIMEPOINT}_sccoda_mcmc_diagnostics.csv", index=False)

az.plot_trace(idata, var_names=["alpha", "beta"])
plt.tight_layout()
plt.savefig(OUT / f"{STUDY}_{TIMEPOINT}_sccoda_mcmc_trace.pdf")
plt.close()

print(diagnostics)
