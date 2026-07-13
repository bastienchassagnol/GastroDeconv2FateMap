"""Harmonise Suppinger and Luque cell-type labels with Cell Ontology (CL).

Workflow
--------
1. Read unmapped H5AD files exported from R via ``anndataR`` (see
   ``scripts/01_02_add_standardised_cell_ontologies.R``).
2. Inspect study-specific labels against CL release ``2023-04-20`` using
   ``bionty.base`` (modern Bionty read-only API; no LaminDB instance required).
3. Apply a **biologist-reviewed crosswalk** that maps each original label to:
   - a broad CL parent term (``cell_type_ontology_term_id``),
   - a neutral cross-study label (``cell_type_harmonized``),
   - a standardised state label (``cell_state_standardized``),
   - a mapping relation describing how literal the match is.
4. Validate assigned CL names and IDs, write crosswalk tables, and save
   harmonised H5AD files for re-import into R as Seurat objects.

Ontology policy
---------------
- **Cell Ontology (CL)** is used for developmental cell types and states.
- **Cell Line Ontology (CLO)** is reserved for experimental lines (H9, H1, …);
  columns are initialised but left empty because source metadata do not encode
  line identity at single-cell level.
- Original author labels are always preserved in ``cell_type_original``.
- Composite, state-level, or artefact labels are mapped to **broad CL parents**
  rather than forced exact matches; ``cell_type_mapping_relation`` records this.

Crosswalk rationale (biologist-guided)
--------------------------------------
Study annotations mix exact CL concepts (e.g. *Epiblast*), germ-layer states
(*Caudal mesoderm* → mesodermal cell), composite slashes (*Epiblast/primitive
streak*), and artefacts (*Zscan4+ Artefact*). The crosswalk keeps NMP and
similar progenitor identities in harmonised columns while assigning defensible
CL parents. Luque ``unknown`` cells remain unresolved (no false CL assignment).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import unicodedata

import anndata as ad
from bionty.base import CellType
import pandas as pd

# =============================================================================
# 0. Paths and ontology version pins
# =============================================================================

ROOT = Path(__file__).resolve().parents[1]
INTERMEDIATE_DIR = ROOT / "data" / "intermediate"
IO = INTERMEDIATE_DIR / "ontology_mapping"
IO.mkdir(parents=True, exist_ok=True)
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

CL_VERSION = "2023-04-20"
CLO_VERSION = "2022-03-21"

# Broad CL parents resolved once against the pinned ontology release.
CL_ONTOLOGY_IDS: dict[str, str] = {
    "embryonic_cell": "CL:0002321",
    "epiblast_cell": "CL:0000352",
    "ectodermal_cell": "CL:0000221",
    "mesodermal_cell": "CL:0000222",
    "endodermal_cell": "CL:0000223",
    "pluripotent_stem_cell": "CL:0002248",
    "neural_cell": "CL:0002319",
    "endothelial_cell": "CL:0000115",
    "paraxial_cell": "CL:0011007",
}

# =============================================================================
# 1. Data structures
# =============================================================================


@dataclass(frozen=True)
class MappingSpec:
    """One row of the biologist-reviewed annotation crosswalk."""

    term_key: str | None
    harmonized: str
    cell_state: str
    relation: str
    artefact: bool = False
    note: str = ""


# =============================================================================
# 2. Crosswalk lookup helpers
# =============================================================================


def annotation_key(value: object) -> str:
    """Normalise labels for case- and accent-insensitive crosswalk lookup."""
    if pd.isna(value):
        return ""
    ascii_text = (
        unicodedata.normalize("NFKD", str(value))
        .encode("ascii", "ignore")
        .decode("ascii")
    )
    return " ".join(ascii_text.strip().casefold().split())


def resolve_cl_terms(cell_type_bt: CellType) -> dict[str, SimpleNamespace]:
    """Resolve pinned CL ontology IDs to preferred names via ``bionty.base``."""
    terms: dict[str, SimpleNamespace] = {}
    for key, ontology_id in CL_ONTOLOGY_IDS.items():
        hits = cell_type_bt.search(ontology_id).reset_index()
        row = hits.loc[hits["ontology_id"].eq(ontology_id)].iloc[0]
        if row["ontology_id"] != ontology_id:
            msg = f"Expected {ontology_id} for {key}, got {row['ontology_id']}."
            raise ValueError(msg)
        terms[key] = SimpleNamespace(name=row["name"], ontology_id=ontology_id)
    return terms


# =============================================================================
# 3. Biologist-reviewed crosswalks
# =============================================================================
# Keys must match ``annotation_key(original_label)``.
# Each entry records: CL parent, harmonised label, state label, and relation.

SUPPINGER_CROSSWALK: dict[str, MappingSpec] = {
    annotation_key("Anterior primitive streak/Def. endoderm"): MappingSpec(
        "embryonic_cell",
        "primitive streak / endoderm",
        "anterior primitive streak / definitive endoderm",
        "composite_to_parent",
        note="Mixed label; mapped to broad embryonic parent.",
    ),
    annotation_key("Caudal epiblast"): MappingSpec(
        "epiblast_cell", "epiblast", "caudal epiblast", "manual_exact_concept"
    ),
    annotation_key("Caudal epiblast/primitive streak"): MappingSpec(
        "epiblast_cell",
        "primitive streak / caudal epiblast",
        "caudal epiblast / primitive streak",
        "composite_to_parent",
    ),
    annotation_key("Caudal mesoderm"): MappingSpec(
        "mesodermal_cell", "mesoderm", "caudal mesoderm", "state_to_parent"
    ),
    annotation_key("Cd63+ ectoderm-like artefact"): MappingSpec(
        "ectodermal_cell",
        "artefact",
        "CD63-positive ectoderm-like artefact",
        "artefact_to_parent",
        artefact=True,
    ),
    annotation_key("Ectopic pluripotency"): MappingSpec(
        "pluripotent_stem_cell",
        "pluripotent",
        "ectopic pluripotency",
        "state_to_parent",
    ),
    annotation_key("Epiblast"): MappingSpec(
        "epiblast_cell", "epiblast", "epiblast", "manual_exact_concept"
    ),
    annotation_key("Epiblast/primitive streak"): MappingSpec(
        "epiblast_cell",
        "primitive streak / caudal epiblast",
        "epiblast / primitive streak",
        "composite_to_parent",
    ),
    annotation_key("Exiting naïve pluripotency"): MappingSpec(
        "pluripotent_stem_cell",
        "pluripotent",
        "exiting naive pluripotency",
        "state_to_parent",
    ),
    annotation_key("Gut"): MappingSpec(
        "endodermal_cell",
        "endoderm",
        "gut-like",
        "state_to_parent",
        note="Tissue/domain label; broad endodermal parent used.",
    ),
    annotation_key("Hemogenic endothelium"): MappingSpec(
        "endothelial_cell",
        "endothelial",
        "hemogenic endothelium",
        "state_to_parent",
    ),
    annotation_key("Naïve pluripotency"): MappingSpec(
        "pluripotent_stem_cell",
        "pluripotent",
        "naive pluripotency",
        "state_to_parent",
    ),
    annotation_key("Neuromesodermal progenitors"): MappingSpec(
        "embryonic_cell",
        "neuromesodermal progenitor",
        "neuromesodermal progenitor",
        "state_to_broad_parent",
        note="NMP identity kept in harmonised columns; no exact CL NMP term.",
    ),
    annotation_key("Paraxial mesoderm"): MappingSpec(
        "paraxial_cell",
        "somitic mesoderm",
        "paraxial mesoderm",
        "synonym_or_near_exact",
    ),
    annotation_key("Pre-somitic mesoderm"): MappingSpec(
        "paraxial_cell",
        "somitic mesoderm",
        "presomitic mesoderm",
        "state_to_parent",
    ),
    annotation_key("Primitive streak"): MappingSpec(
        "embryonic_cell",
        "primitive streak / caudal epiblast",
        "primitive streak",
        "anatomical_state_to_parent",
    ),
    annotation_key("Somite"): MappingSpec(
        "paraxial_cell", "somitic mesoderm", "somite", "anatomical_state_to_parent"
    ),
    annotation_key("Somite differentiation front"): MappingSpec(
        "paraxial_cell",
        "somitic mesoderm",
        "somite differentiation front",
        "state_to_parent",
    ),
    annotation_key("Zscan4+ Artefact"): MappingSpec(
        "pluripotent_stem_cell",
        "artefact",
        "ZSCAN4-positive artefact",
        "artefact_to_parent",
        artefact=True,
    ),
}

LUQUE_CROSSWALK: dict[str, MappingSpec] = {
    annotation_key("caudal epiblast-caudal neuroectoderm"): MappingSpec(
        "embryonic_cell",
        "epiblast / neural",
        "caudal epiblast / caudal neuroectoderm",
        "composite_to_parent",
    ),
    annotation_key("caudal epiblast-primitive streak"): MappingSpec(
        "epiblast_cell",
        "primitive streak / caudal epiblast",
        "caudal epiblast / primitive streak",
        "composite_to_parent",
    ),
    annotation_key("formative"): MappingSpec(
        "pluripotent_stem_cell",
        "pluripotent",
        "formative pluripotency",
        "state_to_parent",
    ),
    annotation_key("meso-biased"): MappingSpec(
        "mesodermal_cell", "mesoderm", "mesoderm-biased", "biased_state_to_parent"
    ),
    annotation_key("naive"): MappingSpec(
        "pluripotent_stem_cell",
        "pluripotent",
        "naive pluripotency",
        "state_to_parent",
    ),
    annotation_key("neural"): MappingSpec(
        "neural_cell", "neural", "neural", "manual_exact_concept"
    ),
    annotation_key("neural-biased"): MappingSpec(
        "neural_cell", "neural", "neural-biased", "biased_state_to_parent"
    ),
    annotation_key("neuromesodermal progenitors"): MappingSpec(
        "embryonic_cell",
        "neuromesodermal progenitor",
        "neuromesodermal progenitor",
        "state_to_broad_parent",
    ),
    annotation_key("pluripotent"): MappingSpec(
        "pluripotent_stem_cell", "pluripotent", "pluripotent", "state_to_parent"
    ),
    annotation_key("primed"): MappingSpec(
        "pluripotent_stem_cell",
        "pluripotent",
        "primed pluripotency",
        "state_to_parent",
    ),
    annotation_key("primitive streak-caudal epiblast"): MappingSpec(
        "epiblast_cell",
        "primitive streak / caudal epiblast",
        "primitive streak / caudal epiblast",
        "composite_to_parent",
    ),
    annotation_key("somitic"): MappingSpec(
        "paraxial_cell", "somitic mesoderm", "somitic", "state_to_parent"
    ),
    annotation_key("unknown"): MappingSpec(
        None,
        "unknown",
        "unknown",
        "unresolved",
        note="No biological identity encoded by source label.",
    ),
}

DEFAULT_UNRESOLVED = MappingSpec(
    None,
    "unknown",
    "unknown",
    "unresolved_not_in_crosswalk",
    note="Label absent from reviewed crosswalk.",
)

# =============================================================================
# 4. Inspection, mapping, validation, and reporting
# =============================================================================


def inspect_labels(
    cell_type_bt: CellType,
    labels: list[str],
    dataset: str,
) -> pd.DataFrame:
    """Run CL ``inspect`` on unique original labels and persist the report."""
    print(f"\n{'=' * 72}\n{dataset}: inspecting {len(labels)} labels\n{'=' * 72}")
    result = cell_type_bt.inspect(labels, field=cell_type_bt.name, return_df=True)
    print(result)
    slug = dataset.lower().replace(" ", "_")
    result.to_csv(IO / f"{slug}_initial_CL_inspection.csv")
    return result


def apply_crosswalk(
    adata: ad.AnnData,
    crosswalk: dict[str, MappingSpec],
    cl_terms: dict[str, SimpleNamespace],
) -> None:
    """Add biological-guided ontology to pre-existing ontology terms."""
    specs = [
        crosswalk.get(annotation_key(v), DEFAULT_UNRESOLVED)
        for v in adata.obs["cell_type_original"].astype("string")
    ]
    records = [cl_terms.get(s.term_key) if s.term_key else None for s in specs]

    adata.obs["cell_type_ontology_name"] = pd.array(
        [r.name if r else pd.NA for r in records], dtype="string"
    )
    adata.obs["cell_type_ontology_term_id"] = pd.array(
        [r.ontology_id if r else pd.NA for r in records], dtype="string"
    )
    adata.obs["cell_type_harmonized"] = pd.array(
        [s.harmonized for s in specs], dtype="string"
    )
    adata.obs["cell_state_standardized"] = pd.array(
        [s.cell_state for s in specs], dtype="string"
    )
    adata.obs["cell_type_mapping_relation"] = pd.array(
        [s.relation for s in specs], dtype="string"
    )
    adata.obs["cell_type_mapping_note"] = pd.array(
        [s.note for s in specs], dtype="string"
    )
    adata.obs["annotation_artefact"] = [s.artefact for s in specs]

    # CLO placeholders — populate only when experimental line names are supplied.
    for col, default in (
        ("cell_line_original", pd.NA),
        ("cell_line_ontology_name", pd.NA),
        ("cell_line_ontology_term_id", pd.NA),
    ):
        if col not in adata.obs:
            adata.obs[col] = pd.array([default] * adata.n_obs, dtype="string")


def validate_cl(cell_type_bt: CellType, adata: ad.AnnData, dataset: str) -> None:
    """Confirm assigned CL names and IDs validate against the pinned release."""
    names = sorted(
        adata.obs["cell_type_ontology_name"].dropna().astype(str).unique().tolist()
    )
    ids = sorted(
        adata.obs["cell_type_ontology_term_id"].dropna().astype(str).unique().tolist()
    )
    print(f"\n{dataset}: final CL name validation")
    print(cell_type_bt.inspect(names, field=cell_type_bt.name, return_df=True))
    print(f"{dataset}: final CL ID validation")
    print(cell_type_bt.inspect(ids, field=cell_type_bt.ontology_id, return_df=True))


def mapping_report(adata: ad.AnnData) -> pd.DataFrame:
    """Summarise cell counts per mapping combination."""
    cols = [
        "dataset",
        "cell_type_original",
        "cell_type_ontology_name",
        "cell_type_ontology_term_id",
        "cell_type_harmonized",
        "cell_state_standardized",
        "cell_type_mapping_relation",
        "cell_type_mapping_note",
        "annotation_artefact",
    ]
    return (
        adata.obs[cols]
        .value_counts(dropna=False)
        .rename("n_cells")
        .reset_index()
        .sort_values(["cell_type_harmonized", "cell_type_original"], kind="stable")
    )


# =============================================================================
# 5. Per-dataset orchestration and pipeline diagram
# =============================================================================


def export_cell_annotations(adata: ad.AnnData, path: Path) -> None:
    """Write per-cell ontology columns to CSV (one row per cell barcode)."""
    cols = [
        "dataset",
        "cell_type_original",
        "timepoint_original",
        "cell_type_ontology_name",
        "cell_type_ontology_term_id",
        "cell_type_harmonized",
        "cell_state_standardized",
        "cell_type_mapping_relation",
        "cell_type_mapping_note",
        "annotation_artefact",
        "cell_line_original",
        "cell_line_ontology_name",
        "cell_line_ontology_term_id",
    ]
    present = [col for col in cols if col in adata.obs.columns]
    df = adata.obs[present].copy()
    df.index.name = "cell_barcode"
    df.to_csv(path)


def harmonise_dataset(
    path: Path,
    out_path: Path,
    crosswalk: dict[str, MappingSpec],
    cell_type_bt: CellType,
    cl_terms: dict[str, SimpleNamespace],
    dataset: str,
) -> ad.AnnData:
    """End-to-end harmonisation for one exported H5AD file."""
    adata = ad.read_h5ad(path)
    labels = sorted(
        adata.obs["cell_type_original"].dropna().astype(str).str.strip().unique()
    )
    labels = [x for x in labels if x]
    inspect_labels(cell_type_bt, labels, dataset)
    apply_crosswalk(adata, crosswalk, cl_terms)
    validate_cl(cell_type_bt, adata, dataset)
    mapping_report(adata).to_csv(
        IO / f"{dataset.lower().replace(' ', '_')}_annotation_crosswalk.csv",
        index=False,
    )
    adata.uns["annotation_harmonization"] = {
        "workflow": "inspect -> reviewed crosswalk -> CL lookup -> validate",
        "cell_type_ontology": "Cell Ontology",
        "cell_type_ontology_version": CL_VERSION,
        "cell_line_ontology": "Cell Line Ontology",
        "cell_line_ontology_version": CLO_VERSION,
        "cell_line_mapping_status": "not populated (no source line names)",
    }
    adata.write_h5ad(out_path)
    return adata


# =============================================================================
# 6. Entry point
# =============================================================================


def main() -> None:
    """Harmonise Suppinger and Luque exports; write annotated H5AD files."""
    cell_type_bt = CellType(source="cl", version=CL_VERSION)
    cl_terms = resolve_cl_terms(cell_type_bt)

    suppinger = harmonise_dataset(
        IO / "suppinger_unmapped.h5ad",
        IO / "suppinger_harmonized.h5ad",
        SUPPINGER_CROSSWALK,
        cell_type_bt,
        cl_terms,
        "Suppinger 2026",
    )
    luque = harmonise_dataset(
        IO / "luque_unmapped.h5ad",
        IO / "luque_harmonized.h5ad",
        LUQUE_CROSSWALK,
        cell_type_bt,
        cl_terms,
        "Luque 2024",
    )

    combined = pd.concat(
        [mapping_report(suppinger), mapping_report(luque)], ignore_index=True
    )
    combined.to_csv(IO / "combined_annotation_crosswalk.csv", index=False)

    export_cell_annotations(
        suppinger,
        INTERMEDIATE_DIR / "suppinger_cell_type_annotations.csv",
    )
    export_cell_annotations(
        luque,
        INTERMEDIATE_DIR / "luque_cell_type_annotations.csv",
    )

    print("\nHarmonised H5AD files written to", IO)
    print("Per-cell annotation CSVs written to", INTERMEDIATE_DIR)


if __name__ == "__main__":
    main()


# =============================================================================
# 7. Helper reference — role and rationale
# =============================================================================
#
# MappingSpec
#   Frozen record for one crosswalk row. Separates the CL parent key from the
#   human-readable harmonised/state labels and documents mapping quality via
#   ``relation`` and optional ``note``. Keeps biology decisions explicit.
#
# annotation_key
#   Normalises author strings before dictionary lookup so accents, spacing, and
#   case differences (e.g. *Naïve* vs *naive*) do not break the crosswalk.
#
# resolve_cl_terms
#   Resolves pinned CL IDs through ``bionty.base`` once at startup. Guarantees
#   ontology names/IDs come from the same CL release used for validation.
#
# SUPPINGER_CROSSWALK / LUQUE_CROSSWALK
#   Biologist-reviewed maps from study-specific labels to broad CL parents.
#   Composite, state, and artefact labels are not forced into false exact CL
#   matches; instead they receive defensible parents plus explicit relations.
#
# inspect_labels
#   Baseline QC step: shows which author labels fail direct CL name validation
#   before crosswalk correction. Output CSV supports audit trails.
#
# apply_crosswalk
#   Cell-level application of reviewed mappings. Adds ontology/harmonised columns
#   while leaving ``cell_type_original`` untouched for traceability.
#
# validate_cl
#   Post-mapping QC: confirms every assigned CL name and ID exists in the
#   pinned ontology release. Catches typos or version drift early.
#
# mapping_report
#   Aggregates cell counts per mapping combination for quick biological review
#   and cross-study comparison tables.
#
# harmonise_dataset
#   Orchestrates the full per-dataset pipeline (inspect → map → validate → export)
#   and stores provenance in ``adata.uns``.
#
# export_cell_annotations
#   Writes one CSV row per cell with harmonised ontology metadata for downstream
#   R/Python workflows outside AnnData.
#
# main
#   Runs Suppinger and Luque harmonisation, writes combined crosswalk CSV, and
#   saves the pipeline diagram alongside other ontology-mapping artefacts.
