# GastroDeconv2FateMap

Using a combination of bulk transcriptome, single-cell, imaging, and
deconvolution algorithms, predict the phenotypic endpoint of gastruloids
from early stages

## Run the package with VSCode custom keybindings

| Shortcut       | Action                                         |
|----------------|------------------------------------------------|
| `Ctrl+Shift+B` | **R: devtools::load_all** (default build task) |
| `Ctrl+Shift+D` | **R: devtools::document**                      |
| `Ctrl+Shift+T` | **R: devtools::test**                          |
| `Ctrl+Shift+C` | **R: devtools::check**                         |

## Data management with DVC

Large datasets in this project are tracked with [DVC](https://dvc.org/)
and stored on a shared remote. Git holds small `.dvc` pointer files; the
actual data files live on the remote at
`/mnt/DATA_11TB/projects/dtoo_project/dvc-remotes/GastroDeconv2FateMap`.

Tracked locations:

- `data/raw/` — raw inputs (e.g. GEO downloads, Seurat objects)
- `data/intermediate/` — processed objects produced by analysis scripts

### Prerequisites

Install DVC with SSH support:

``` bash
uv install "dvc[ssh]"
```

Clone the repository and activate the R environment (`renv`) as usual.
DVC must be available in your shell before running the commands below.

### 1. Configure access to the DVC remote

The project default remote is `mmg_cluster`. The shared configuration
(in `.dvc/config`) points to:

``` text
ssh://mmg-sb-05:/mnt/DATA_11TB/projects/dtoo_project/dvc-remotes/GastroDeconv2FateMap
```

Each user must set up **machine-local** settings in `.dvc/config.local`
(this file is not committed to Git).

#### Option A — SSH remote (cluster or remote workstation)

**Step 1 — Create a dedicated SSH key** (do not copy private keys
between machines):

``` bash
mkdir -p ~/.ssh && chmod 700 ~/.ssh
ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519_mmg -C "yourname@your-host" -N ""
chmod 600 ~/.ssh/id_ed25519_mmg
chmod 644 ~/.ssh/id_ed25519_mmg.pub
```

**Step 2 — Authorise the public key** on the storage host. Ask a project
maintainer to add your `~/.ssh/id_ed25519_mmg.pub`, or run (if you have
access):

``` bash
ssh-copy-id -i ~/.ssh/id_ed25519_mmg.pub <user>@<storage-host>
```

Replace `<storage-host>` with a hostname that resolves from your
machine. If `mmg-sb-05` fails with *Temporary failure in name
resolution*, use the resolvable alias instead (e.g. `S-T-MMG-SB-05`).

**Step 3 — Test SSH access:**

``` bash
ssh -i ~/.ssh/id_ed25519_mmg <user>@<storage-host> \
  "ls /mnt/DATA_11TB/projects/dtoo_project/dvc-remotes/GastroDeconv2FateMap"
```

**Step 4 — Configure DVC locally** (creates or updates
`.dvc/config.local`):

``` bash
dvc remote modify --local mmg_cluster \
  url ssh://<user>@<storage-host>:/mnt/DATA_11TB/projects/dtoo_project/dvc-remotes/GastroDeconv2FateMap

dvc remote modify --local mmg_cluster \
  keyfile ~/.ssh/id_ed25519_mmg
```

**Step 5 — Verify the remote:**

``` bash
dvc remote list
dvc status
```

#### Option B — Local path remote (storage already mounted)

If the remote directory is already visible on your filesystem (same NFS
mount or local path), skip SSH entirely:

``` bash
dvc remote modify --local mmg_cluster \
  url /mnt/DATA_11TB/projects/dtoo_project/dvc-remotes/GastroDeconv2FateMap
```

### 2. Download tracked datasets (`dvc pull`)

After cloning or pulling Git changes that include new `.dvc` files:

``` bash
git pull
dvc pull
```

Pull only a specific file or folder:

``` bash
dvc pull data/raw/GSE229513_gastruloidsobject.rds.dvc
dvc pull data/raw/
dvc pull data/intermediate/
```

Check which files are missing locally before pulling:

``` bash
dvc status
```

### 3. Add and publish new datasets (`dvc add` + `dvc push`)

Place new files under `data/raw/` or `data/intermediate/`, then track
them with DVC. The data file itself stays out of Git; only the `.dvc`
pointer is committed.

**Step 1 — Track one or more files:**

``` bash
dvc add data/raw/my_new_dataset.rds
```

> DVC creates a `.dvc` stub next to each file and updates
> `data/raw/.gitignore` or `data/intermediate/.gitignore` automatically.

**Step 2 — Upload data to the remote:**

``` bash
dvc push
```

**Step 3 — Commit the DVC pointers to Git:**

``` bash
git add data/raw/my_new_dataset.rds.dvc data/raw/.gitignore
git commit -m "Track my_new_dataset with DVC"
git push
```

**Other collaborators can then run `git pull` followed by `dvc pull` to
obtain the new files.**

### Troubleshooting

| Symptom | Likely cause | Fix |
|----|----|----|
| `Permission denied (publickey)` | SSH key not authorised on storage host | Add your public key; set `keyfile` in `.dvc/config.local` |
| `Connection timed out` | Network cannot reach storage host | Use Option B if storage is mounted locally, or contact IT for VPN/routing |
