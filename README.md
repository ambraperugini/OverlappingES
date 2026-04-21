# OverlappingES

Repository accompanying the paper:

> *Perugini, A., Calignano, G., & Pastore, M. (2026)*  
> **Not normal: a simulation study comparing effect sizes for skewed psychological data**  
> Frontiers in Psychology

---

## Overview

This repository contains all materials used in the study, including:

- Simulation code  
- Data generation procedures  
- Figures used in the manuscript  
- Revision documents  
- Final published paper  

The project investigates how different **effect size indices behave under non-normal conditions**, focusing on skewness and variance heterogeneity.

---

## Repository Structure


OverlappingES/
│
├── Rcodes/ # Core R scripts for simulations and analysis
├── data/ # Data used/generated in the study
├── knitr/ # Reproducible report components (R Markdown / Quarto)
├── Submission figures/ # Figures included in the manuscript
│
├── fpsyg-17-1771029.pdf # Published paper
├── Revision_001.pdf # Revision document
├── Revision_001.Qmd # Source for revision
├── Comments_to_Reviewers.qmd
│
├── overlapping.bib # Bibliography
├── apa.csl # Citation style
├── license.txt # CC BY 4.0 license
├── README.md
└── .gitignore


---

## Scientific Content

The repository implements a **simulation-based comparison** of four effect size indices:

- **Cohen’s d**  
- **Common Language Effect Size (CLES)**  
- **Parametric overlap (ηp)**  
- **Non-parametric overlap (η)**  

These are evaluated under controlled violations of:

- Normality  
- Variance homogeneity  
- Symmetry  

---

## Simulation Design

- Reference distribution: Normal (mean = 0, variance = 1)  
- Comparison distribution: Skew-normal  

Manipulated parameters:

- Mean difference (δ): 0, 2  
- Variance (σ): 1, 5  
- Skewness (α): 0, 10  
- Sample sizes: 10 → 1000  
- Replications: 2000 per condition  

---

## Key Results

- **Cohen’s d**
  - Most stable and reliable across conditions  
  - Robust to violations  

- **CLES and ηp**
  - Highly correlated with d  
  - Poorer inferential performance under non-normality  

- **η (overlap index)**
  - Captures full distribution differences  
  - Robust to skewness and heteroscedasticity  
  - Complements mean-based indices  

Core result:

> Effect sizes are not interchangeable, even when strongly correlated.

---

## Usage

### Requirements

- R (recommended version ≥ 4.0)  
- Required packages (see scripts in `Rcodes/`)

### Run analysis

Execute scripts inside:


Rcodes/


### Reproduce manuscript

Use files inside:


knitr/


---

## Reproducibility

The repository is structured to support:

- Full simulation replication  
- Figure regeneration  
- Transparent methodological workflow  

---

## License

This project is licensed under:

**Creative Commons Attribution 4.0 International (CC BY 4.0)**

See `license.txt` for full terms.

---

## Authors

- Ambra Perugini  
- Giulia Calignano  
- Massimiliano Pastore  

University of Padua

---
## Notes
This is research-oriented code:
- Not optimized for production  
- Designed for methodological transparency  
- Focused on reproducibility and statistical validation
- Designed for methodological transparency  
- Focused on reproducibility and statistical validation
