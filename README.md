# Project Overview

The folder contains the essential code of the data analysis associated with the Bermejo Ruiz, M. et al, "Micro-Cm: restrictase-free microbiome-wide chromosome conformation profiling", 2026 manuscript.

---

## Folder structure:

- `code/`:
	- Defined_community_plas_chrom_contacts.ipynb - Code for figures 3b and 3c.
	- Obtention_of_contact_profiles.ipynb - Code for figure 2.
	- Scailing_curves.ipynb - Code for figure 3a.
  - Stool_sample_secondary_analysis.ipynb - Code for figure 4a and supplementary figure 4.
  - fragment_genomes.py - split of reference genomes into contigs
  - run_hicexplorer.sh - commands to run HicExplorer and construct contact matrices (like the ones in figure 1 and supplementary figure 2)

- `download_data.sh`: code for downloading from Zenodo the input files for the notebooks.

---

## Instructions:

For installing the repository and downloading the input data please run:
```bash
git clone https://github.com/leylabmpi/microc_paper_code
cd microc_paper_code
./download_data.sh
```
The installation is expected to be complete in a few minutes.

The R notebooks (.ipynb files) are provided as prerendered versions (viewable via a Web-browser at Github) - a part of the visualizations and other results presented in the manuscript. 

Then install the following non-Conda -based packages into it:
```bash
R
devtools::install_github("leylabmpi/LeyLabRMisc")
```

Open a notebook in VS Code, select a R Jupyter kernel and run the notebook.
