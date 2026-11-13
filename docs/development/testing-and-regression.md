# Testing and regression checks

> Maintainer reference for the v0.3.0 release.

Run checks from the repository root in an isolated environment:

```bash
python -m pip install -e .
python -m pip install 'pytest>=8,<9' 'ruff>=0.14,<1'
python -m ruff check src tests
python -m compileall -q src
python -m pytest -q
```

The published model-only example is also a smoke test:

```bash
mkdir -p results
moselecta \
  --inputs tests/data/alternative_splice_smoke_a.gff3 tests/data/alternative_splice_smoke_b.gff3 \
  --source-name predictor_a predictor_b \
  --out-prefix results/example
test -s results/example.final_annotation.gff3
```

Documentation verification also checks local Markdown links, documented option
parsing/defaults, output keys/suffixes, DOCX/PDF rendering, distribution
contents, clean-source extraction, and installation from the wheel. Scientific
invariance is checked by comparing `src/` and `tests/` byte-for-byte with the
accepted baseline and, where available, the established representative fixture
without normalising outputs.

For a documentation-only change, any source or test diff is a release blocker
unless separately reviewed and authorized. Generated caches, environments,
build directories, and internal verification reports do not belong in public
source distributions.

The editable guide lives in `docs/user-guide/`. After a material edit,
regenerate its PDF, render every page to images, inspect the complete document,
and verify the table of contents and PDF bookmarks before replacement.

