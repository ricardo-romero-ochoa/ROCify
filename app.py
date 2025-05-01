#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import streamlit as st
from data_processing import run_data_pipeline
from roc_analysis import run_roc_analysis
import pandas as pd

st.title("Automated ROC Analysis for TCGA Genes")

# User Inputs
tcga_project = st.selectbox("Select TCGA Project", ["TCGA-LIHC", "TCGA-BRCA", "TCGA-GBM", "..."])
genes_input = st.text_input("Enter genes (comma-separated)", "CDK1,AURKB,CCNA2")

if st.button("Run Analysis"):
    genes = [g.strip() for g in genes_input.split(",")]

    with st.spinner("Processing Data..."):
        run_data_pipeline(tcga_project, genes)
        run_roc_analysis(genes)

    st.success("Analysis Complete!")

    st.image("outputs/roc_curves.png", caption="ROC Curves")

    auc_df = pd.read_csv("outputs/auc_results.csv")
    st.subheader("AUC Results")
    st.write(auc_df)
    st.download_button("Download AUC Results", auc_df.to_csv(index=False), "auc_results.csv")

    cutoff_df = pd.read_csv("outputs/optimal_cutoffs.csv")
    st.subheader("Optimal Cutoffs")
    st.write(cutoff_df)
    st.download_button("Download Optimal Cutoffs", cutoff_df.to_csv(index=False), "optimal_cutoffs.csv")

