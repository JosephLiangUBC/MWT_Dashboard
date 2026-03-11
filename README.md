# MWT_Dashboard

Data dashboard for Multi-Worm Tracker (MWT) experiments from the Rankin Lab.

## Overview

This repository contains a Streamlit dashboard for exploring consolidated *C. elegans* MWT datasets, along with the notebooks and utilities used to process data before upload to the backing database.

The current dashboard:

- loads experiment data from a PostgreSQL backend
- supports dataset-level filtering across screens
- provides gene-level and allele-level phenotype exploration
- supports multi-gene and multi-allele comparison views
- includes post-stimulus arousal (PSA) summaries
- includes a clustering page for exploratory grouping of genes

Youtube overview:

[![MWT Dashboard Introduction Video](http://img.youtube.com/vi/xbQ0CRUCnZs/0.jpg)](https://youtu.be/xbQ0CRUCnZs "MWT Dashboard Introduction Video")

## Repository Layout

- `MWT_dashboard.py`: current Streamlit entrypoint
- `app_pages/`: dashboard pages
- `utils/`: data loading, preprocessing, auth, and helper functions
- `config.py`: shared plotting configuration and page metadata
- `Data_processing/`: data-processing workflow notebooks
- `MWT_whitepaper.md`: project and workflow documentation
- `legacy/`: older dashboard and analysis assets retained for reference

## Running the Dashboard Locally

The current app entrypoint is `MWT_dashboard.py`.

```bash
streamlit run MWT_dashboard.py
```

The app expects Streamlit secrets for authentication and database access. Create `.streamlit/secrets.toml` with the values used by `utils/auth.py` and `utils/data_loader.py`:

```toml
password = "..."
psql_user = "..."
psql_passwword = "..."
```

Notes:

- the dashboard connects to the `mwtdata` PostgreSQL database on AWS RDS
- without valid secrets, the app will not pass the login gate or load data
- the default local Streamlit port is `8501`

## Installing Dependencies

### Using pip

1. Clone the repository.
2. Create and activate a virtual environment.

```bash
python -m venv rankinlab
source rankinlab/bin/activate     # macOS/Linux
rankinlab\Scripts\activate.bat    # Windows
```

3. Install the required dependencies.

```bash
pip install -r requirements.txt
```

4. If dependencies change, refresh the lockfile.

```bash
pip freeze > requirements.txt
```

### Using conda

1. Clone the repository.
2. Create and activate the conda environment.

```bash
conda env create -f environment.yml
conda activate rankinlab
```

3. If the environment definition changes, update it.

```bash
conda env update --file environment.yml --prune
```

## Docker

The repository also includes a `Dockerfile` that starts the current Streamlit app on port `8502`.

```bash
docker build -t mwt-dashboard .
docker run --rm -p 8502:8502 mwt-dashboard
```

If you use Docker, make sure the container has access to the required Streamlit secrets.

## Dashboard Pages

- `Home - Getting Started`: project background, data description, and usage guidance
- `Data at a Glance`: phenotype ranking and whole-dataset heatmap views
- `Gene-specific Data`: phenotype profile, ranking, response curves, and exports for one gene
- `Allele-specific Data`: phenotype profile, ranking, response curves, and exports for one allele
- `Custom Gene Selection`: compare multiple genes side by side
- `Custom Allele Selection`: compare multiple alleles side by side
- `Post Stimulus Data`: PSA summaries across selected genes and metrics
- `Gene Clustering`: PCA + KMeans + t-SNE exploratory clustering
- `Citations`: references for the dashboard and source datasets

## Data and Processing Notes

- The dashboard is designed around the 10 s ISI, 30-tap habituation protocol described in the app.
- Heatmaps and phenotypic profiles are derived from normalized t-stat style summaries, with N2 controls used as the reference within screen.
- Gene and allele profile tables are generated dynamically in `utils/data_loader.py` from backend tables rather than being stored as static dashboard assets.
- The notebooks in this repository remain part of the workflow for preparing and validating upload-ready data.

## New Changes Since the Last README Update

The last README-focused update landed on **May 26, 2025**. Since then, the dashboard has changed substantially:

- the active app entrypoint is now `MWT_dashboard.py`; the old `MWT_dashboard_copy.py` is now under `legacy/`
- password protection was added and then reimplemented through Streamlit session-state auth
- the dashboard expanded from the original single-gene/allele views to include:
  - `Custom Gene Selection`
  - `Custom Allele Selection`
  - `Post Stimulus Data`
  - `Gene Clustering`
- PSA metrics were added and then overhauled across the pipeline and dashboard
- heatmap/profile generation moved toward backend-driven loading and in-app transformation
- download support was expanded with static plot exports and CSV exports across pages
- Docker support was added for running the dashboard in a container
- the help and whitepaper documentation were expanded
- 2026 updates introduced additional backend and heatmap formatting changes, including support for updated t-score / corrected p-value table formats

In short, the repository is now a broader dashboard-and-pipeline workspace rather than just a basic Streamlit viewer.

## Contact

For bugs, issues, or feature requests, contact Joseph Liang at `joseph.liang@psych.ubc.ca`.
