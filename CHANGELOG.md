# Changelog

All notable changes to this project will be documented in this file.

## [v1.9] - 2025-06-04

New release of nf-core/yascp, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- Unified pipeline to use a **single container** with consistent software versions.
- Container built without reliance on paid Conda channels.
- Cleaned and annotated all config files.
- Standardized input: Doublet and Celltype modules now expect `.h5ad` only.
- Moved `celltype` assignment earlier in the pipeline for modular reusability.

### `Fixed`
- Inconsistent input expectations across modules (e.g. `.mtx` vs `.h5ad`).
- Minor bugs affecting results folder consistency.

### `Dependencies`
- Nextflow >=24.10.4
- Singularity >=3.11.4
- Removed dependency on user-managed Conda environments.

### `Deprecated`
- Dropped support for `.mtx` inputs in `--doublets` and `--celltypes` modes.


## [v1.8] - 2025-02-01

### `Added`
- Hashtag-based deconvolution module.
- Improved folder structure for results.
- Output of final Seurat-compatible `.rds` files.

### `Fixed`
- Path resolution issues and sample-level subdirectory outputs.

### `Dependencies`
- Same as v1.7

### `Deprecated`
- N/A


## [v1.7] - 2024-11-01

Used for all *Cardinal Freeze* datasets.

### `Added`
- Doublet detection support for any `AnnData` object.
- Improved logging and diagnostics in doublet mode.

### `Fixed`
- Memory handling for large h5ad files.

### `Dependencies`
- Same as v1.6.1

### `Deprecated`
- N/A


## [v1.6.1] - 2024-08-01

### `Added`
- Improved logic for conditional module execution.

### `Fixed`
- Bugs in restart behavior, directory caching, and CLI wrapper consistency.

### `Dependencies`
- Same as v1.6

### `Deprecated`
- N/A


## [v1.6] - 2024-07-01

### `Added`
- Fully modular execution across ambient RNA, celltype, clustering, and doublets.
- Consistent h5ad/matrix file compatibility.

### `Fixed`
- Input format validation and improved test mode.

### `Dependencies`
- Adds `bcftools 1.0`, `irods 4.2.7`, `baton` support.

### `Deprecated`
- N/A


## [v1.5] - 2024-05-01

### `Added`
- Five additional doublet detection methods.
- Entry points: `--doublets`, `--cellbender`, user-supplied CellTypist and Azimuth references.

### `Fixed`
- Module-specific logs and improved error reporting.
- VDJ and Citeseq compatibility improvements.

### `Dependencies`
- Introduced additional module loads for iRODS and BCFtools.

### `Deprecated`
- N/A


## [v1.4] - 2024-03-01

### `Added`
- Cell Ranger 7 compatibility.
- Option to skip integration during reclustering.
- GPU-free test mode (disables CellBender).
- Better Citeseq filtering.

### `Fixed`
- Logging cleanup and edge-case handling in Citeseq mode.

### `Dependencies`
- Same as v1.3

### `Deprecated`
- N/A


## [v1.3] - 2023-12-01

### `Added`
- `--cluster` mode (`--JUST_RECLUSTER`) added for subclustering and integration-free workflows.

### `Fixed`
- Output labelling and subfolder creation.

### `Dependencies`
- Same as v1.2

### `Deprecated`
- N/A


## [v1.2] - 2023-10-01

### `Added`
- Improved documentation and user messages for `--JUST_CELLTYPES`.

### `Fixed`
- Bugfixes in test mode and entry point handling.

### `Dependencies`
- Same as v1.1

### `Deprecated`
- N/A


## [v1.1] - 2023-08-01

### `Added`
- New entry point: `--celltypes` to run only cell type assignment.

### `Fixed`
- Minor issues in path resolution and version messaging.

### `Dependencies`
- Same as v1.0

### `Deprecated`
- N/A


## [v1.0] - 2023-06-01

Initial release.

### `Added`
- Core pipeline for processing 10X scRNA-seq data: ambient RNA removal, donor deconvolution, celltype assignment, genotype concordance.
- Multiple run modes: `test`, `fetch`, `sample_input`, `clean`, full `input.nf`.

### `Fixed`
- N/A

### `Dependencies`
- Singularity 3.11.4
- Nextflow 22.04.4

### `Deprecated`
- N/A
