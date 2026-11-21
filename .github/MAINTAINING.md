# Maintaining BindFlow

Notes for whoever administers the repository. Contributors want `CONTRIBUTING.md` instead.

## Codacy

Most of the analysis is configured in the repository, but some of it only exists in the
Codacy interface, so it is written down here.

### Set in the Codacy interface

Under **Code patterns** of the repository:

- **Ruff**: enabled, with *"Use Ruff configuration file found in the repository"* turned
  on. That makes Codacy read `[tool.ruff]` from `pyproject.toml` of the branch being
  analysed instead of its own list of patterns. The file has to be pushed for it to be
  found: on a branch without it, Codacy falls back to the patterns that were set before the
  option was turned on.
- **Prospector**: its docstring patterns turned off, or the tool disabled altogether, now
  that Ruff covers the same ground.

### Why Prospector is turned down

Codacy has no pydocstyle tool of its own: those checks come from Prospector, which bundles
pydocstyle, and Codacy enables **D212** and **D213** together. The two state the opposite
of each other, "the summary of a multi-line docstring goes on the first line" against "on
the second line", so with both enabled every multi-line docstring is reported whatever it
does. There are 81 of them in `src/bindflow` and the count stays at 81 however they are
written. The package follows D212, the style PEP 257 describes.

Ruff is the only Python linter whose configuration file Codacy reads, so none of this can
be moved into the repository.

### What `.codacy.yml` can do

`.codacy.yml` at the root excludes paths from the analysis, and that is all it does here:
it cannot enable or disable tools or rules. What it excludes is prose, `docs/` and every
`.md` and `.rst`, so that a page of text is not reported as code, plus the deprecated
folder. The code stays in, the dashboard and the contributed scripts included.

## Cutting a release

1. `docs/source/CHANGELOG.md`: `## [unreleased]` becomes `## [X.Y.Z] - <date>`, with the
   compare link at the bottom of the file next to the others. The tags carry a `v` since
   v3.4.1, and the links have to match, or they 404.
2. Merge the release branch into `main` without squashing, so the history survives.
3. Tag `main` after the merge, not the branch: versioningit derives the version from the
   tag, and a tag on the branch would sit behind the merge commit.

   ```bash
   git checkout main && git pull
   git tag -a vX.Y.Z -m "bindflow X.Y.Z"
   ```

4. Build and check before publishing anything, because the GitHub release triggers the
   upload on its own:

   ```bash
   python -m build
   twine check dist/*
   ```

5. Push the tag, then create the GitHub release from it. That is what runs
   `.github/workflows/pypi-publish.yml`.
