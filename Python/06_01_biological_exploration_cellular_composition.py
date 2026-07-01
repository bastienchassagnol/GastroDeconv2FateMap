"""Smoke test for pertpy in the project uv virtual environment."""

# %%
import importlib.metadata

import pertpy as pt
a = 20; b = 30; c = a + b

# %%

def main() -> None:
    print(f"pertpy {pt.__version__}")
    print(f"scanpy {importlib.metadata.version('scanpy')}")

    adata = pt.dt.dialogue_example()
    print(
        "Loaded pertpy dialogue example: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )
    print("pertpy installation OK")


if __name__ == "__main__":
    main()
