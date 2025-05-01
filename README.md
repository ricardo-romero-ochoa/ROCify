# ROCify: Automated TCGA Gene Expression Analysis

**ROCify** is an automated pipeline designed to simplify the process of querying gene expression data from The Cancer Genome Atlas (TCGA). The application efficiently performs ROC (Receiver Operating Characteristic) analysis, computes Area Under the Curve (AUC) values, and identifies optimal expression cut-offs to facilitate biomarker discovery.

## Features

- **Data Acquisition**: Seamlessly download and process RNA-seq data from any TCGA project.
- **ROC Curve Generation**: Automatically generate high-quality ROC curves for user-specified genes.
- **AUC Calculation**: Accurately compute AUC values to evaluate biomarker potential.
- **Optimal Cut-off Determination**: Quickly determine optimal gene expression thresholds using Youden's Index.

## Getting Started

Clone this repository to your local machine:

```bash
git clone https://github.com/your-username/rocify.git
cd rocify
```

### Requirements

Ensure you have R and Python installed. Install required packages:

```bash
pip install -r requirements.txt
```
```bash
Rscript install_packages.R
```

### Usage

Run the Streamlit app locally:

```bash
streamlit run app.py
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

## License

Distributed under the GPLv3 License.

