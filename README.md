# Epigenetic Architectures of Chronic Pain: A Multi-Pipeline Methylation Analysis of Fibromyalgia

![R](https://img.shields.io/badge/R-4.3.0-blue.svg)
![Pipeline](https://img.shields.io/badge/Pipeline-ChAMP%20%7C%20Minfi-orange.svg)
![Status](https://img.shields.io/badge/Status-Completed-success.svg)

## 📌 Executive Summary
Fibromyalgia is a debilitating chronic pain disorder with an elusive etiology. This project challenges the "functional" nature of the disease by uncovering **robust epigenetic alterations** in the structural and neurological machinery of patients.

Using a rigorous **Multi-Pipeline Consensus Approach (ChAMP + Minfi)** on Illumina 450k array data, I identified **43 high-confidence consensus genes** that withstood strict cross-validation. The findings propose a novel **"Dual-Hit Pathology"**:
1.  **Structural Fragility:** Epigenetic silencing of *TNXB* (Ehlers-Danlos gene) and *AGRN*.
2.  **Synaptic Rewiring:** Maladaptive plasticity driven by *NRXN1* and the *PCDH* cluster.

---

## 🚀 The "Consensus" Methodology
Unlike standard analyses that rely on a single algorithm, I implemented a **comparative benchmarking strategy** to eliminate false positives.

| Pipeline Step | Tool Used | Purpose | Outcome |
| :--- | :--- | :--- | :--- |
| **Exploration** | **ChAMP** | Broad signal detection (Linear Regression) | Identified ~18,000 DMPs (High Sensitivity) |
| **Validation** | **Minfi** | Regional clustering (Bump Hunting) | Identified 120 DMRs (High Specificity) |
| **Consensus** | **Custom R** | **Intersection Analysis** | **43 Robust Biomarkers** (Genes found by BOTH methods) |

> *By filtering the broad ChAMP signal through the strict Minfi regional filter, I isolated the "true biological signal" from technical noise.*

---

## 🔬 Key Biological Discoveries

### 1. The "Synaptic Velcro" Defect
Functional network analysis (STRING-db) of the 43 consensus genes revealed a distinct cluster of **Synaptic Adhesion Molecules**.
* **The Genes:** `NRXN1` (Neurexin-1), `CTNNA2` (Catenin), `PCDHA7` (Protocadherin).
* **The Insight:** These proteins physically hold the pre- and post-synaptic neurons together. Their epigenetic disruption suggests **"Maladaptive Synaptic Plasticity"**—the neurons are physically rewiring themselves to maintain a hypersensitive pain state.

### 2. The "Ehlers-Danlos" Connection
I identified significant hypermethylation in **Tenascin XB (`TNXB`)**.
* **Clinical Relevance:** Mutations in *TNXB* cause Ehlers-Danlos Syndrome (EDS).
* **The Insight:** This provides a molecular link explaining the high comorbidity between Fibromyalgia and Joint Hypermobility. The patient's extracellular matrix (ECM) may be structurally compromised.

### 3. Epigenetic Feedback Loops
The consensus list included **`HDAC4`** (Histone Deacetylase 4) and **`YTHDC1`** (RNA Methylation Reader).
* **The Insight:** The disease signature includes the very machinery that controls gene expression, suggesting a self-perpetuating **"Epigenetic Trap"** that prevents recovery.

---

## 📊 Visualizations

### A. The Consensus Strategy
The Venn diagram demonstrates the strict filtering process. Only genes validated by both the Probe-Level (ChAMP) and Region-Level (Minfi) analysis were selected.
![Venn Diagram](results/figures/Venn_Consensus_Genes.png)

### B. Protein-Protein Interaction (PPI) Network
The 43 consensus genes form tight functional modules regulating **Synaptic Structure** (Red) and **RNA Processing** (Green).
![PPI Network](results/figures/PPI_Network_Consensus.png)

### C. Target Validation: TNXB
Methylation levels of *TNXB*, showing clear separation between Fibromyalgia and Control samples.
![TNXB Boxplot](results/figures/TNXB_Boxplot.png)

---

## 🛠️ Tech Stack & Skills Demonstrated
* **Genomic Data Science:** `ChAMP`, `minfi`, `GenomicRanges`, `limma`.
* **Statistical Rigor:** False Discovery Rate (FDR) correction, SVD Batch Effect Correction, Consensus Calling.
* **Network Biology:** PPI analysis (STRING), Functional Enrichment (DAVID/KEGG).
* **Visualization:** `ggplot2`, `VennDiagram`, `ComplexHeatmap`.

---

## 📬 Contact
**Sivananth M.**
*Bioinformatics Analyst & Epigenetics Enthusiast*
* [LinkedIn](https://www.linkedin.com/in/sivananth-m)
* [Email](mailto:sivananthm2002@gmail.com)

> *"I believe that the future of precision medicine lies in distinguishing signal from noise. This project is my demonstration of that belief."*
