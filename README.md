# metabolomics-training

[![Tests](https://github.com/benson1231/metabolomics-training/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/benson1231/metabolomics-training/actions/workflows/test.yml)
[![Docker Engine](https://img.shields.io/badge/Docker-27.5.1-blue?logo=docker)](https://www.docker.com/get-started/)


This repository provides a minimal, reproducible workflow to practice **metabolomics data analysis** using **MetaboAnalystR 4.0** inside a Docker container.

The goal is to offer a clean, environment‑isolated setup so users can focus on LC-MS preprocessing, statistical analysis, enrichment/pathway analysis, and functional interpretations.

---

## 🚀 Run MetaboAnalystR with Docker

Make sure Docker and Docker Compose are installed. Then run:

```bash
docker compose run --rm metaboanalystr4
```

This launches an isolated R environment with all necessary packages preinstalled.

---

## 📂 Run Example Workflows (inside the container)

All analysis scripts are located in the `scripts/` directory. Execute them inside the container as follows:

```bash
Rscript scripts/introduction.R

Rscript scripts/read_data.R

Rscript scripts/lcms.R

Rscript scripts/functional_analysis.R

Rscript scripts/statistical.R

Rscript scripts/enrichment.R
```

All results will be automatically saved to the `results/` directory.

---

## 📘 Official MetaboAnalystR 4.0 Tutorials

Below is a curated list of helpful references:

* [MetaboAnalystR 4.0 Official Tutorial](https://www.metaboanalyst.ca/docs/RTutorial.xhtml)
* [Introduction](https://www.metaboanalyst.ca/resources/vignettes/Introductions.html)
* [LC-MS/MS Raw Spectra Processing](https://www.metaboanalyst.ca/resources/vignettes/LCMSMS_Raw_Spectral_Processing.html)
* [Functional Analysis of Global Metabolomics](https://www.metaboanalyst.ca/resources/vignettes/Functional_Analysis_global_metabolomics.html)
* [Statistical Analysis (one-factor)](https://www.metaboanalyst.ca/resources/vignettes/Statistical_Analysis_Module.html)
* [Enrichment Analysis of Targeted Metabolomics](https://www.metaboanalyst.ca/resources/vignettes/Enrichment_Analysis.html)
* [Pathway Analysis of Targeted Metabolomics](https://www.metaboanalyst.ca/resources/vignettes/Pathway_Analysis.html)
* [Biomarker Analysis](https://drive.google.com/file/d/1DjLDE9IGvU_rjfdIDdlUtS7FGahfcsov/view)
* [Statistical Analysis (Metadata table)]()
* [Joint-Pathway Analysis](https://drive.google.com/file/d/1HVoXNX98CZLcpr7MkpVZP5DueXcnZeIj/view)
* [Functional Meta-Analysis]()
* [Network Analysis](https://drive.google.com/file/d/1fZ364APP8pqemFg0oD0kIhAzH9Zm9p-0/view)
* [Power Analysis](https://drive.google.com/file/d/1UE7592V1vvpaeeyCfVPD_692Coxww238/view)
* [Meta-Analysis](https://drive.google.com/file/d/1nCH980icadGcQJwtQ8yHDKcuM32V7BLG/view)

