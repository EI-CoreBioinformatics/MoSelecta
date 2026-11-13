# Installing MoSelecta

MoSelecta is pure Python and has no runtime dependencies. Install it with `pip`, or run it from the Docker image without installing Python at all.

| Method | Command | You need |
|---|---|---|
| pip | `pip install moselecta` | Python 3.10 to 3.15 |
| Docker | `docker pull gemygk/moselecta` | Docker (or Apptainer/Singularity) |

## Requirements

- for the pip route: CPython `3.10`, `3.11`, `3.12`, `3.13`, `3.14` or `3.15`, and `pip` in an isolated environment
- for the Docker route: Docker, or Apptainer/Singularity; Python is not needed
- no required runtime Python dependencies beyond the standard library
- write access to the chosen environment and output directories

Installing a released version with `pip` does not need a build. A source build through `pip` or `python -m build` normally creates an isolated build environment and downloads the build backend named in `pyproject.toml` (`uv_build`), so it may require package-index access unless the backend is already cached or provided locally. Offline source builds are not supported.

## Install with pip (in a virtual environment)

A virtual environment keeps MoSelecta separate from your other Python packages.

On Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install moselecta
```

On Windows PowerShell:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
pip install moselecta
```

To install a specific release, or to upgrade an existing installation:

```bash
pip install moselecta==0.3.0
pip install --upgrade moselecta
```

## Run with Docker

The image [`gemygk/moselecta`](https://hub.docker.com/r/gemygk/moselecta) on Docker Hub contains MoSelecta and Python, so nothing else needs installing. It is built for `linux/amd64` and `linux/arm64`.

```bash
docker pull gemygk/moselecta:0.3.0
docker run --rm gemygk/moselecta:0.3.0 moselecta --version
```

The `latest` tag always points to the newest release. Use a version tag, as above, when a pipeline must stay reproducible.

On HPC systems that provide Apptainer or Singularity instead of Docker, use the same image directly:

```bash
apptainer exec docker://gemygk/moselecta:0.3.0 moselecta --version
```

## Install a wheel

Install from a wheel file you already have:

```bash
python -m pip install /path/to/moselecta-0.3.0-py3-none-any.whl
```

The wheel is platform-independent because MoSelecta is pure Python.

## Install a source distribution

```bash
python -m pip install /path/to/moselecta-0.3.0.tar.gz
```

Build isolation resolves the compatible `uv_build` backend unless it is already available from the configured environment/cache.

## Install from a source checkout

From the root of a cloned or unpacked copy of the repository:

```bash
python -m pip install .
```

For editable development work:

```bash
python -m pip install -e .
```

Editable installation is for development and test work, not the final release artifact.

## Build wheel and source distribution

Install a PEP 517 build frontend, then build both artifacts:

```bash
python -m pip install "build>=1,<2"
python -m build
```

The outputs are written to `dist/` by default. The frontend creates an isolated build environment and uses the backend constraint in `pyproject.toml`. A build performed with another frontend is acceptable only if it honours the same PEP 517 configuration.

## Check the installation

```bash
moselecta --version
moselecta --help
moselecta --full-help
moselecta --expert-help
moselecta-compare --help
```

`moselecta --full-help` and `moselecta --expert-help` list the less common options. The companion command `moselecta-compare --help` is installed alongside.

## Try a small example

This confirms that MoSelecta runs end to end; it is not a real analysis. The two small input files are in the source repository, not in the installed package, so clone the repository first:

```bash
git clone https://github.com/EI-CoreBioinformatics/moselecta.git
cd moselecta
mkdir -p results
moselecta \
  --inputs tests/data/alternative_splice_smoke_a.gff3 tests/data/alternative_splice_smoke_b.gff3 \
  --source-name predictor_a predictor_b \
  --out-prefix results/example
```

`--out-prefix` is a filename prefix, not a directory-creation option. Its parent (`results` above) must already exist and be writable.

A successful run ends with `Run complete` and writes `results/example.final_annotation.gff3`, which contains one gene with one transcript.

With Docker, run the same command from the cloned directory with `docker run` in front:

```bash
mkdir -p results
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data gemygk/moselecta:0.3.0 \
  moselecta \
  --inputs tests/data/alternative_splice_smoke_a.gff3 tests/data/alternative_splice_smoke_b.gff3 \
  --source-name predictor_a predictor_b \
  --out-prefix results/example
```

The container works in `/data`. Mount the directory that holds your inputs there and give paths relative to it.

`-u "$(id -u):$(id -g)"` makes the output files belong to you instead of to root.

## Run a realistic example

For a complete demonstration, the [Bombus 4 Mb example](examples/README.md) provides seven gene-prediction sources, genomic sequence, protein and transcript evidence, and validated MoSelecta outputs. It illustrates a realistic integration workflow, including optional UTR augmentation, and provides expected results for checking a new installation.

## Before your first run

- **Create the output directory.** `--out-prefix results/run1` writes files such as `results/run1.final_annotation.gff3`. MoSelecta does not create `results/`; make it first and check that it is writable. Use a new prefix for each run so that old outputs are not mistaken for new ones.
- **Allow for disk space with the larger output profiles.** The default `standard` profile is right for routine use. `--output-profile full` and `debug` write a large set of audit files.
- **Give `--genome-fasta` an uncompressed FASTA.** MoSelecta reads the index `<fasta>.fai` next to it, creates the index if it is missing, and rebuilds it if it is out of date, so the FASTA directory must be writable in those cases. On shared storage, build the index once in advance and keep it. Sequence names in the annotations must match the FASTA sequence names.

## Development verification

For a contributor checkout, install the package and the development tools, then run the checks:

```bash
python -m pip install -e .
python -m pip install "pytest>=8,<9" "ruff>=0.14,<1"
python -m pytest
python -m ruff check src tests
python -m compileall -q src tests
```

The release-preparation process additionally inspects wheel/sdist contents, installs the wheel in a fresh environment, executes published examples, checks documentation links/options/schemas, and renders the user guide.

## Troubleshooting

**`pip` rejects the Python version.** MoSelecta supports Python 3.10 to 3.15. Check with `python --version`, or use the Docker image.

**`moselecta: command not found`.** Activate the environment you installed into, then check with `pip show moselecta`.

**Output files cannot be opened.** Create the parent directory of `--out-prefix` and make sure it is writable. With Docker, the directory must be inside the folder mounted at `/data`.

**Docker writes files owned by root.** Add `-u "$(id -u):$(id -g)"` to `docker run`.

**The FASTA index cannot be created.** Make the FASTA directory writable, or put a valid `.fai` next to the FASTA.

**A source build cannot find its build backend.** Check network access to the package index. Installing a released version with `pip install moselecta` does not need a build.

## Related documentation

- [README](README.md)
- [Configuration](docs/configuration.md)
- [Outputs](docs/outputs.md)
- [Bombus 4 Mb example](examples/README.md)
- [MoSelecta user guide](docs/user-guide/MoSelecta-v0.3.0-user-guide.pdf)
