# Spatial Transcriptomics Analysis — Tutorial Notes

A consolidated reference guide covering key concepts, tools, and workflows from four spatial transcriptomics tutorials using **Scanpy** and **Squidpy**.

---

## Table of Contents

1. [Overview](#overview)
2. [Tutorial 1 — Scanpy: Basic Spatial Analysis of Visium Data](#tutorial-1--scanpy-basic-spatial-analysis-of-visium-data)
3. [Tutorial 2 — Squidpy: Analyzing Visium Fluorescence Data](#tutorial-2--squidpy-analyzing-visium-fluorescence-data)
4. [Tutorial 3 — Squidpy: Analyzing Visium H&E Data](#tutorial-3--squidpy-analyzing-visium-he-data)
5. [Tutorial 4 — Squidpy: Analyzing Xenium Data](#tutorial-4--squidpy-analyzing-xenium-data)
6. [Key Concepts Learned](#key-concepts-learned)
7. [Common Tools & Functions Reference](#common-tools--functions-reference)
8. [Resources & Links](#resources--links)

---

## Overview

Spatial transcriptomics allows us to measure gene expression while preserving the physical location of cells within a tissue. These tutorials explore how to load, preprocess, visualize, and analyze spatial data using two core Python libraries:

- **[Scanpy](https://scanpy.readthedocs.io/)** — A scalable toolkit for analyzing single-cell gene expression data, with dedicated spatial visualization functions.
- **[Squidpy](https://squidpy.readthedocs.io/)** — Built on top of Scanpy/AnnData, Squidpy extends the ecosystem with image analysis, spatial graph construction, and spatial statistics tools.

All data is stored in the **AnnData** (`.h5ad`) format, which holds count matrices, metadata, spatial coordinates, and images in a unified structure.

---

## Tutorial 1 — Scanpy: Basic Spatial Analysis of Visium Data

**Source:** https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html

### Dataset

A 10x Genomics **Visium** spatial transcriptomics dataset of the **human lymph node**, publicly available from the 10x Genomics website. The dataset was loaded using:

```python
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
```

> For custom Visium data, use `squidpy.read.visium()` or `sc.read_visium()` to import it.

### What We Learned

#### Data Loading & Structure
- Visium data contains **UMI counts**, **high-resolution tissue images**, and **spatial spot coordinates** packaged together in an AnnData object.
- Each Visium spot captures transcripts from approximately 5–20 cells, making it spatially resolved but not single-cell resolution.
- The AnnData object for Visium includes `adata.obsm["spatial"]` for coordinates and `adata.uns["spatial"]` for image data.

#### Quality Control (QC)
- Standard QC metrics are computed using `sc.pp.calculate_qc_metrics()`.
- Mitochondrial gene fraction is calculated and used as a QC indicator:

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
```

- Spots with very low total counts or very few expressed genes are filtered out to remove poor-quality observations.

#### Spatial Visualization
- `sc.pl.spatial()` overlays gene expression or cluster labels directly onto the H&E tissue image.
- Key parameters:
  - `img_key` — specifies which image resolution to use (e.g., `"hires"`)
  - `crop_coord` — crops the image to a region of interest
  - `alpha_img` — controls image transparency
  - `size` — acts as a **scaling factor** for spot sizes (different from non-spatial plots)

```python
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
```

#### Normalization & Clustering
- Counts are normalized using `sc.pp.normalize_total()`, followed by log-transformation (`sc.pp.log1p()`).
- Highly variable genes (HVGs) are identified with `sc.pp.highly_variable_genes()`.
- Dimensionality reduction: **PCA** → **UMAP**.
- Clustering with **Leiden algorithm** (`sc.tl.leiden()`).
- Clustering was performed in **gene expression space**, and results were visualized in both UMAP and spatial coordinates to gain insight into **tissue organization**.

#### Spatial Gene Expression Analysis
- By mapping clusters onto the tissue image, distinct spatial domains such as germinal centers, mantle zones, and T-cell zones can be identified in the lymph node.
- This highlights how **spatial context** reveals patterns that UMAP alone cannot show.

---

## Tutorial 2 — Squidpy: Analyzing Visium Fluorescence Data

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html

### Dataset

A pre-processed Visium slide of a **coronal section of the mouse brain**, provided as an AnnData object with pre-annotated clusters, paired with a **fluorescence tissue image** in `squidpy.im.ImageContainer` format.

```python
img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()
```

The fluorescence image contains three channels:
- **DAPI** — specific to DNA (nucleus staining)
- **anti-NEUN** — specific to neurons
- **anti-GFAP** — specific to glial cells

### What We Learned

#### ImageContainer
- Squidpy introduces `squidpy.im.ImageContainer` as the standard data structure for tissue images. It supports multi-channel fluorescence images and enables cropping, processing, and feature extraction per Visium spot.
- Image channels can be visualized individually with `img.show(channelwise=True)`.

#### Image Preprocessing
- Before segmentation, it is best practice to smooth the image using `sq.im.process()` to reduce noise.

#### Image Segmentation
- `sq.im.segment()` performs nucleus segmentation using the **DAPI channel**.
- Segmentation masks can then be used to count nuclei per Visium spot, providing a **cell density estimate**.

#### Image Feature Extraction
- `sq.im.calculate_image_features()` extracts quantitative features from the tissue image for each Visium spot, creating an `obs × features` matrix stored in `adata.obsm`.
- Features extracted in this tutorial:
  - **Summary features** — statistical summaries (mean, std) of pixel intensity per spot
  - **Histogram features** — pixel intensity distributions
  - **Segmentation features** — nucleus count, mean nucleus area
  - **Texture features** — GLCM-based texture descriptors
- Features can be computed at **multiple crop sizes and scales** to capture multi-scale context.

#### Feature-Based Clustering
- Extracted image features were used to compute an independent **morphology-based cluster annotation** using the Leiden algorithm.
- This cluster annotation was compared to the gene expression–based clustering to find **correlated and complementary information**.
- Key insight: image morphology and gene expression carry **overlapping** information (e.g., distinct cell types have distinct morphology) as well as **complementary** information (e.g., cell density from segmentation is not directly in the expression matrix).

#### Spatial Visualization
- `sq.pl.spatial_scatter()` is Squidpy's primary function for plotting clusters or gene expression overlaid on the tissue image.

---

## Tutorial 3 — Squidpy: Analyzing Visium H&E Data

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html

### Dataset

A pre-processed Visium slide of a **coronal section of the mouse brain** (H&E stained), provided with pre-annotated clusters. This tutorial focuses on spatial statistics and graph analysis in addition to image features.

```python
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()
```

### What We Learned

#### Multi-Scale Image Features
- Summary features were extracted at two different scales (`scale=1.0` and `scale=2.0`) to capture both local and broader tissue context.
- Features from multiple scales were concatenated into a single `adata.obsm["features"]` matrix for downstream clustering.

```python
for scale in [1.0, 2.0]:
    sq.im.calculate_image_features(adata, img.compute(), features="summary",
                                   key_added=f"features_summary_scale{scale}", scale=scale)
```

#### Morphology-Based Clustering
- Image features alone can be used to cluster spots by their **tissue morphology**, independent of gene expression.
- This demonstrates that H&E images encode biologically meaningful spatial structure.

#### Spatial Graph Construction
- A **spatial neighborhood graph** is built using `sq.gr.spatial_neighbors()`, which connects each Visium spot to its neighbors based on physical proximity.
- This graph is the foundation for all spatial statistics analyses.

#### Spatial Statistics

**Co-occurrence Score:**
- `sq.gr.co_occurrence()` calculates how often pairs of cell types (clusters) co-occur at various spatial distances.
- This answers: *do certain cell types tend to appear near each other?*

**Ripley's Statistics:**
- `sq.gr.ripley()` computes Ripley's L or K function to assess spatial clustering or dispersion of a cell type.
- A value above the theoretical expectation indicates **spatial clustering**.

**Moran's I (Spatial Autocorrelation):**
- `sq.gr.spatial_autocorr()` computes Moran's I to identify genes whose expression is spatially patterned.
- High Moran's I = gene expression is **spatially autocorrelated** (nearby spots have similar expression).
- This is used to identify spatially variable genes (SVGs).

**Neighborhood Enrichment:**
- `sq.gr.nhood_enrichment()` tests whether certain cluster pairs are **over- or under-represented** as direct neighbors on the spatial graph compared to random expectation.

#### Visualization
- `sq.pl.co_occurrence()`, `sq.pl.ripley()`, `sq.pl.spatial_scatter()` provide publication-ready spatial plots.

---

## Tutorial 4 — Squidpy: Analyzing Xenium Data

**Source:** https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html

### Dataset

A **Xenium** (10x Genomics) dataset of **human lung cancer** tissue, representing sub-cellular resolution spatial transcriptomics. Data was loaded using `spatialdata-io`:

```python
from spatialdata_io import xenium
sdata = xenium(xenium_path)
sdata.write(zarr_path)
```

The dataset contains **161,000 cells** and **480 genes**.

### What We Learned

#### SpatialData Framework
- Xenium data is handled via the **SpatialData** object (`spatialdata` library), which organizes:
  - **Images** — morphology focus image (multi-resolution)
  - **Labels** — cell and nucleus segmentation masks
  - **Points** — individual transcript locations (~40 million)
  - **Shapes** — cell and nucleus boundary polygons
  - **Tables** — AnnData object with count matrix and metadata
- Data is stored in the efficient **Zarr** format for scalable read/write.

#### Data Ingestion
- The `spatialdata-io` library provides a `xenium()` reader that automatically parses and converts Xenium outputs to the SpatialData format.
- The AnnData table is extracted from `sdata.tables["table"]` for standard Scanpy/Squidpy analysis.

#### Xenium-Specific Metadata
- Each cell has rich metadata including: `transcript_counts`, `control_probe_counts`, `cell_area`, `nucleus_area`, `z_level`, and `nucleus_count`.
- This metadata is critical for quality filtering and normalization decisions.

#### Quality Control
- QC metrics specific to Xenium include monitoring `control_probe_counts` (should be near zero) alongside standard `total_counts` and expressed gene filters.
- `sc.pp.calculate_qc_metrics()` is applied as in standard Scanpy workflows.

#### Spatial Coordinates
- Cell coordinates are stored in `adata.obsm["spatial"]` (automatically set by `spatialdata-io`).
- For manually prepared data, this must be set explicitly.

#### UMAP & Spatial Visualization
- Standard Scanpy UMAP workflow is applied, and results are visualized both in UMAP space and in spatial coordinates using `sq.pl.spatial_scatter()`.

#### Spatial Statistics
- The same spatial statistics tools from Tutorial 3 apply here (neighborhood enrichment, co-occurrence, Ripley's, Moran's I), but at **single-cell resolution** rather than the Visium spot level.

#### Interactive Visualization with napari-spatialdata
- `napari-spatialdata` is introduced as an interactive viewer that allows exploration of spatial data including transcript-level point clouds, segmentation masks, and cluster annotations in a GUI environment.

---

## Key Concepts Learned

### AnnData Structure for Spatial Data
| Slot | Content |
|---|---|
| `adata.X` | Count matrix (cells × genes) |
| `adata.obs` | Cell/spot metadata (QC metrics, cluster labels) |
| `adata.var` | Gene metadata (mitochondrial flag, etc.) |
| `adata.obsm["spatial"]` | 2D spatial coordinates |
| `adata.uns["spatial"]` | Tissue image and scale factors |
| `adata.obsm["features"]` | Extracted image features |

### Spatial Resolution Comparison
| Platform | Resolution | Cells per spot | Tutorial |
|---|---|---|---|
| Visium | ~55 µm spots | 5–20 cells | Tutorial 1, 2, 3 |
| Xenium | Sub-cellular | Single cell | Tutorial 4 |

### Analysis Workflow (General)
1. Load data (AnnData / SpatialData)
2. QC and filtering
3. Normalization and HVG selection
4. Dimensionality reduction (PCA → UMAP)
5. Clustering (Leiden)
6. Spatial visualization (overlay on tissue image)
7. Image feature extraction (Squidpy)
8. Spatial graph construction
9. Spatial statistics (co-occurrence, Moran's I, neighborhood enrichment)

---

## Common Tools & Functions Reference

### Scanpy
| Function | Purpose |
|---|---|
| `sc.datasets.visium_sge()` | Download and load a Visium dataset |
| `sc.read_visium()` | Load custom Visium data |
| `sc.pp.calculate_qc_metrics()` | Compute QC statistics |
| `sc.pp.normalize_total()` | Library size normalization |
| `sc.pp.highly_variable_genes()` | Identify highly variable genes |
| `sc.tl.leiden()` | Leiden clustering |
| `sc.pl.spatial()` | Spatial overlay visualization |
| `sc.pl.umap()` | UMAP visualization |

### Squidpy
| Function | Purpose |
|---|---|
| `sq.im.ImageContainer` | Multi-channel tissue image container |
| `sq.im.process()` | Image preprocessing (e.g., smoothing) |
| `sq.im.segment()` | Nucleus segmentation |
| `sq.im.calculate_image_features()` | Extract image features per spot |
| `sq.gr.spatial_neighbors()` | Build spatial neighborhood graph |
| `sq.gr.nhood_enrichment()` | Neighborhood enrichment test |
| `sq.gr.co_occurrence()` | Spatial co-occurrence score |
| `sq.gr.ripley()` | Ripley's spatial statistics |
| `sq.gr.spatial_autocorr()` | Moran's I / spatial autocorrelation |
| `sq.pl.spatial_scatter()` | Spatial scatter visualization |

### SpatialData / Xenium
| Function | Purpose |
|---|---|
| `spatialdata_io.xenium()` | Read Xenium output directory |
| `sdata.write(zarr_path)` | Save SpatialData to Zarr format |
| `sd.read_zarr()` | Load SpatialData from Zarr |
| `sdata.tables["table"]` | Access AnnData count table |

---

## Resources & Links

| Resource | URL |
|---|---|
| Tutorial 1 — Scanpy Basic Spatial Analysis | https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html |
| Tutorial 2 — Squidpy Visium Fluorescence | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html |
| Tutorial 3 — Squidpy Visium H&E | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html |
| Tutorial 4 — Squidpy Xenium | https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html |
| Scanpy Documentation | https://scanpy.readthedocs.io/ |
| Squidpy Documentation | https://squidpy.readthedocs.io/ |
| AnnData Documentation | https://anndata.readthedocs.io/ |
| SpatialData Documentation | https://spatialdata.scverse.org/ |
| 10x Genomics Datasets | https://support.10xgenomics.com/spatial-gene-expression/datasets |
| Allen Brain Atlas | https://mouse.brain-map.org/ |

---

*This README was compiled from the four tutorials listed above and summarizes the key concepts, tools, and workflows covered across them.*
