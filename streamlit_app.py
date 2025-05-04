
#!/usr/bin/env python
# coding: utf-8

# ------------------------------------------------------------------------------
# ROCify: Automated ROC Analysis for TCGA Genes
# © 2025 DataBiotica. All rights reserved.
# ------------------------------------------------------------------------------


import streamlit as st
from data_processing import run_data_pipeline
from roc_analysis import run_roc_analysis
from PIL import Image
import pandas as pd
import os

# ✅ Must be first
st.set_page_config(
    page_title="ROCify",
    page_icon="💻",
    layout="wide"
)

subprocess.run(["Rscript", "install_packages.R"], check=True)


# ✅ Branding
logo = Image.open("assets/logo.png")
st.image(logo, width=150)

st.title("💻 ROCify: Automated ROC Analysis for TCGA Genes")
st.caption("ROC curves and optimal cutoffs analyses for biomarker prognosis.")


# TCGA Project Selection with Descriptions
tcga_projects = {
    "ACC – Adrenocortical carcinoma": "TCGA-ACC",
    "BLCA – Bladder Urothelial Carcinoma": "TCGA-BLCA",
    "BRCA – Breast invasive carcinoma": "TCGA-BRCA",
    "CESC – Cervical squamous cell carcinoma and endocervical adenocarcinoma": "TCGA-CESC",
    "CHOL – Cholangio carcinoma": "TCGA-CHOL",
    "COAD – Colon adenocarcinoma": "TCGA-COAD",
    "DLBC – Lymphoid Neoplasm Diffuse Large B-cell Lymphoma": "TCGA-DLBC",
    "ESCA – Esophageal carcinoma": "TCGA-ESCA",
    "GBM – Glioblastoma multiforme": "TCGA-GBM",
    "HNSC – Head and Neck squamous cell carcinoma": "TCGA-HNSC",
    "KICH – Kidney Chromophobe": "TCGA-KICH",
    "KIRC – Kidney renal clear cell carcinoma": "TCGA-KIRC",
    "KIRP – Kidney renal papillary cell carcinoma": "TCGA-KIRP",
    "LAML – Acute Myeloid Leukemia": "TCGA-LAML",
    "LGG – Brain Lower Grade Glioma": "TCGA-LGG",
    "LIHC – Liver hepatocellular carcinoma": "TCGA-LIHC",
    "LUAD – Lung adenocarcinoma": "TCGA-LUAD",
    "LUSC – Lung squamous cell carcinoma": "TCGA-LUSC",
    "MESO – Mesothelioma": "TCGA-MESO",
    "OV – Ovarian serous cystadenocarcinoma": "TCGA-OV",
    "PAAD – Pancreatic adenocarcinoma": "TCGA-PAAD",
    "PCPG – Pheochromocytoma and Paraganglioma": "TCGA-PCPG",
    "PRAD – Prostate adenocarcinoma": "TCGA-PRAD",
    "READ – Rectum adenocarcinoma": "TCGA-READ",
    "SARC – Sarcoma": "TCGA-SARC",
    "SKCM – Skin Cutaneous Melanoma": "TCGA-SKCM",
    "STAD – Stomach adenocarcinoma": "TCGA-STAD",
    "TGCT – Testicular Germ Cell Tumors": "TCGA-TGCT",
    "THCA – Thyroid carcinoma": "TCGA-THCA",
    "THYM – Thymoma": "TCGA-THYM",
    "UCEC – Uterine Corpus Endometrial Carcinoma": "TCGA-UCEC",
    "UCS – Uterine Carcinosarcoma": "TCGA-UCS",
    "UVM – Uveal Melanoma": "TCGA-UVM"
}

# Initialize session state
if "analysis_complete" not in st.session_state:
    st.session_state.analysis_complete = False

if st.button("🔄 Reset and Clear Results"):
    st.session_state.analysis_complete = False
    for f in ["outputs/roc_curves.png", "outputs/auc_results.csv", "outputs/optimal_cutoffs.csv", "expr_data.csv"]:
        if os.path.exists(f):
            os.remove(f)
    st.experimental_rerun()

selected_label = st.selectbox("Choose a cancer type:", list(tcga_projects.keys()))
tcga_project = tcga_projects[selected_label]
genes_input = st.text_input("Enter one or more genes (comma-separated), below are just examples", "CDK1,AURKB,CCNA2")

if st.button("Run Analysis"):
    genes = [g.strip() for g in genes_input.split(",")]
    st.info("⚠️ Downloading and processing TCGA data can take several minutes. Please be patient.")

    with st.status("Step 1 of 2: Downloading and processing TCGA data. This may take a few minutes...", expanded=True):
        st.write("🔄 Contacting GDC API and downloading raw counts...")
        run_data_pipeline(tcga_project, genes)

    with st.status("Step 2 of 2: Performing ROC analysis...", expanded=True):
        st.write("📊 Calculating AUC and generating plots...")
        run_roc_analysis(genes)

    st.session_state.analysis_complete = True
    st.success("Analysis Complete!")

if st.session_state.analysis_complete:
    if os.path.exists("outputs/roc_curves.png"):
        st.image("outputs/roc_curves.png", caption="ROC Curves")
    else:
        st.warning("ROC plot not found.")

    if os.path.exists("outputs/auc_results.csv"):
        auc_df = pd.read_csv("outputs/auc_results.csv")
        st.subheader("AUC Results")
        st.write(auc_df)
        st.download_button("Download AUC Results", auc_df.to_csv(index=False), "auc_results.csv")

    if os.path.exists("outputs/optimal_cutoffs.csv"):
        cutoff_df = pd.read_csv("outputs/optimal_cutoffs.csv")
        st.subheader("Optimal Cutoffs")
        st.write(cutoff_df)
        st.download_button("Download Optimal Cutoffs", cutoff_df.to_csv(index=False), "optimal_cutoffs.csv")


# ----------------------------------------------------------------------
# Footer
# ----------------------------------------------------------------------
st.markdown("---")
st.caption("© 2025 DataBiotica / ROCify. All rights reserved.")
