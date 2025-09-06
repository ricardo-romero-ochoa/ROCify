# ROCify: Automated TCGA Gene Expression Analysis

[![DOI](https://zenodo.org/badge/975904628.svg)](https://doi.org/10.5281/zenodo.15347297)

**ROCify** is an automated pipeline designed to simplify the process of querying gene expression data from The Cancer Genome Atlas (TCGA). The application efficiently performs ROC (Receiver Operating Characteristic) analysis, computes Area Under the Curve (AUC) values, and identifies optimal expression cut-offs to facilitate biomarker discovery.

## Features

- **Data Acquisition**: Seamlessly download and process RNA-seq data from any TCGA project.
- **ROC Curve Generation**: Automatically generate high-quality ROC curves for user-specified genes.
- **AUC Calculation**: Accurately compute AUC values to evaluate biomarker potential.
- **Optimal Cut-off Determination**: Quickly determine optimal gene expression thresholds using Youden's Index.

## Getting Started

Clone this repository to your local machine:

```bash
git clone https://github.com/ricardo-romero-ochoa/rocify.git
cd rocify
```

### Requirements

Ensure you have R and Python installed. For Windows, make sure Rscript is in your PATH

Alternatevely for Windows users, edit this line in data_processing.py and streamlit_app.py

```bash
rscript_path = shutil.which("Rscript") or "C:\\Program Files\\R\\R-4.5.1\\bin\\Rscript.exe"
```
to point to the Rscript executable

**Python packages**:

Run the following command directly from your terminal (for Windows, use Anaconda prompt):

```bash
pip install -r requirements.txt
```

**R packages**:

Run the following command from your terminal (you may need admin privileges):

```bash
Rscript install_packages.R
```
Alternatevely, install Bioconductor packacges manually from the terminal:

```bash
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
and then:

```bash
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "edgeR"))
```

### Usage

Run the Streamlit app locally (for Windows, use Anaconda prompt):

```bash
streamlit run streamlit_app.py
```

Then, open your browser and navigate to `http://localhost:8501`.

## Output

The app provides:
- Downloadable ROC curves plots
- Tabulated AUC values
- Optimal cut-off tables

## Applications

- Biomarker discovery
- Cancer genomics research
- Diagnostic tool development

## Contributing

Contributions are welcome! Submit issues and pull requests to improve this tool.

## Testing
The app has been tested on Ubuntu 24.04.2 and on Windows 11

## License

Distributed under the GPLv3 License.

