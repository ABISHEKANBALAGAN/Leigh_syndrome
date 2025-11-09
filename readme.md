# Leigh Syndrome Proteomics Analysis Pipeline

## 📋 Overview

A comprehensive bioinformatics pipeline for re-analyzing proteomics data from Leigh Syndrome patient-derived neural models. This multi-day workflow transforms raw mass spectrometry data into prioritized therapeutic targets and biomarkers through integrated statistical, network, and functional analysis.

**🔬 Study Context**: Based on "Defective metabolic programming impairs early neuronal morphogenesis in neural cultures and an organoid model of Leigh syndrome" (Inak et al., 2021). This pipeline re-analyzes the proteomics data to identify dysregulated pathways and prioritize biomarkers.

## 🎯 Pipeline Summary

| Day | Analysis | Key Outputs |
|-----|----------|-------------|
| 0 | Data Preparation | Cleaned protein groups, metadata |
| 1 | Preprocessing | Normalized, imputed intensity matrix |
| 2 | Differential Analysis | 1,101 significant proteins |
| 3 | Functional Enrichment | 1,330 enriched pathways |
| 4 | Network Analysis | 154-protein PPI network |
| 5 | Biomarker Prioritization | Ranked therapeutic targets |

## 📁 Project Structure

```
Leigh_syndrome/
├── data/
│   ├── proteinGroups.txt              # Raw MaxQuant output
│   ├── proteinGroups_clean.txt        # Filtered contaminants
│   └── metadata.tsv                   # Sample information
├── scripts/
│   ├── 0_processing/
│   │   └── preprocess.py              # Data extraction
│   ├── 1_preprocessing/
│   │   └── preprocess_filter_normalize.py  # QC & normalization
│   ├── 2_DAanalysis/
│   │   └── DAA.py   # Statistical testing
│   ├── 3_Functional_enrichment_and_Biological_annotation/
│   │   └── FEBA.py   # Pathway analysis
│   ├── 3_Functional_enrichment_and_Biological_annotation/
│   │   └── string_ppi.R  # PPI networks
│   └── 4_target_biomarker/
│       └── biomarker_prioritization.py # Multi-criteria scoring
├──lfq_matrix_intensities_only.csv     # Data Extraction outputs 
├──lfq_matrix_with_annotations.csv     # Data Extraction outputs
├──protein_info.csv                    # Data Extraction outputs
├── preprocessing_results/             # Day 1 outputs
│   ├── 01_log2_transformation.png
│   ├── 02_missing_value_diagnostics.png
│   ├── 03_normalization_effects.png
│   ├── 04_imputation_effects.png
│   ├── 05_correlation_pca.png
│   ├── intensities_imputed.csv
│   ├── preprocessed_data.pkl
│   ├── preprocessed_data_full.csv
│   └── protein_info_filtered.csv
├── differential_results/              # Day 2 outputs
│   ├── differential_results_full.csv
│   ├── heatmap_top_proteins.png
│   ├── ma_plot.png
│   ├── pvalue_distribution.png
│   ├── qc_summary.csv
│   ├── sig_for_enrichment.tsv
│   ├── significant_proteins.csv
│   ├── top100_proteins.csv
│   └── volcano_plot.png
├── v3_enrichment_results/             # Day 3 outputs
│   ├── biological_summary_annotated.csv
│   ├── enrichment_all_dotplot.png
│   ├── enrichment_down_dotplot.png
│   ├── enrichment_up_dotplot.png
│   ├── gprofiler_all_significant.csv
│   ├── gprofiler_down_regulated.csv
│   ├── gprofiler_up_regulated.csv
│   ├── mitochondrial_proteins.csv
│   ├── pathway_comparison.png
│   └── top100_annotated_summary.csv
├── new_1_enrichment_results/          # Day 4 outputs
│   ├── cytoscape_edges.txt
│   ├── cytoscape_nodes.txt
│   ├── string_enrichment_component.csv
│   ├── string_enrichment_function.csv
│   ├── string_enrichment_kegg.csv
│   ├── string_enrichment_pfam.csv
│   ├── string_enrichment_process.csv
│   ├── string_gene_mapping.csv
│   ├── string_input_genes.txt
│   ├── string_interactions.csv
│   ├── string_modules.csv
│   ├── string_network_communities.png
│   ├── string_network_hubs.png
│   ├── string_network_regulation.png
│   ├── string_network_statistics.csv
│   └── string_network_statistics_plot.png
└── biomarker_results/                 # Day 5 outputs
    ├── biomarker_heatmap_top30.png
    ├── biomarker_prioritization_summary.png
    ├── biomarker_ranking_full.csv
    ├── biomarker_report.txt
    ├── biomarker_top50.csv
    ├── biomarker_volcano_prioritized.png
    ├── known_leigh_genes_validated.csv
    └── validation_checklist.txt
```

## 🛠️ Installation & Dependencies

### Python Environment
```bash
# Create conda environment
conda create -n leigh_proteomics python=3.9
conda activate leigh_proteomics

# Install Python packages
pip install pandas numpy matplotlib seaborn scipy statsmodels scikit-learn requests
```

### R Environment
```bash
# Install required R packages
Rscript -e "install.packages('ggplot2', 'dplyr', 'readr', 'igraph', repos='https://cloud.r-project.org')"
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('STRINGdb')"
```

## 🚀 Usage

### Day 0: Data Preparation
```bash
# Filter contaminants and create metadata
grep -v "CON__" data/proteinGroups.txt | grep -v "REV__" > data/proteinGroups_clean.txt

# Create metadata from summary file
awk 'BEGIN { OFS="\t"; print "RawFile","Experiment","Fraction","Group","Condition","Replicate" }
     NR>1 && $1 ~ /_/ && $3 ~ /^[0-9]+$/ && $1 !~ /^Total$/ {
       split($1,a,"_");
       group=a[1];
       rep=a[2];
       cond = (group=="3P" ? "Patient" : "Control");
       print $1, $2, $3, group, cond, rep
     }' data/summary.txt > data/metadata.tsv
```

### Day 1: Data Extraction & Preprocessing
```bash
python scripts/0_processing/preprocess.py
python scripts/1_preprocessing/preprocess_filter_normalize.py
```
**Outputs**: Normalized log2 intensities, QC plots, imputed data matrix

### Day 2: Differential Expression Analysis
```bash
python scripts/2_DAanalysis/DAA.py
```
**Outputs**: 1,101 significant proteins, volcano plots, statistical results

### Day 3: Functional Enrichment
```bash
python scripts/3_Functional_enrichment_and_Biological_annotation/FEBA.py
```
**Outputs**: Pathway enrichments, mitochondrial analysis, biological summaries

### Day 4: Network Analysis
```bash
Rscript scripts/3_Functional_enrichment_and_Biological_annotation/string_ppi.R
```
**Outputs**: PPI network, hub proteins, functional modules, Cytoscape files

### Day 5: Biomarker Prioritization
```bash
python scripts/4_target_biomarker/biomarker_prioritization.py
```
**Outputs**: Ranked biomarkers, multi-criteria scores, validation checklist

## 📊 Key Results

### Statistical Findings
- **1,101 significantly altered proteins** (21.2% of proteome)
- **Balanced regulation**: 599 up, 502 down in controls
- **Strong effects**: Median 2.1× fold change(logFC = 1.049) in significant proteins
- **High confidence**: Empirical Bayes moderation with prior df = 1.00, prior variance = 0.0455
- **Stringent thresholds**: |logFC| > 1.0 (2-fold change), adj.P.Val < 0.05 (FDR correction)

### Biological Insights
- **Nuclear reprogramming**: DNA replication, ribosome biogenesis upregulated
- **Extracellular communication**: Vesicles, secretion pathways downregulated  
- **Mitochondrial vulnerability**: 74 mitochondrial proteins affected
- **Complex-specific effects**: Complex I and IV most impacted

### Network Topology
- **154 connected proteins** with 1,018 high-confidence interactions
- **Hub proteins**: PARP1 (64 connections), CD44 (48), MCM complex members
- **Modular structure**: 9 functional communities (modularity = 0.596)

### Top Biomarker Candidates
1. **NDUFAF5** (Score: 57.2) - Known Leigh gene, Complex I assembly
2. **NDUFA4** (Score: 45.3) - Complex I subunit
3. **MT-CO2** (Score: 43.1) - Mitochondrial-encoded Complex IV
4. **NDUFV3** (Score: 42.4) - Complex I subunit
5. **COX6B1** (Score: 39.9) - Complex IV assembly

## 🔬 Validation

### Cytoscape Network Confirmation
The pipeline generates ready-to-use Cytoscape files:
- `cytoscape_nodes.txt` - Protein attributes and scores
- `cytoscape_edges.txt` - Interaction network

**Import to Cytoscape**:
1. File → Import → Network from File (edges.txt)
2. File → Import → Table from File (nodes.txt)
3. Visualize hub proteins and functional modules

### Experimental Validation Checklist
The pipeline provides a comprehensive wet-lab validation protocol including:
- Western blot targets and antibody recommendations
- qPCR primer design list
- Functional assay protocols (complex activities, ATP production)
- Sample size calculations and statistical considerations
- Budget estimates ($10,000-15,000 for full validation)

## 💡 Interpretation Guide

### For Biologists
- Focus on **Tier 1 & 2 candidates** for immediate validation
- Prioritize **mitochondrial Complex I/IV proteins**
- Consider **known Leigh genes** (NDUFAF5, SURF1 homologs)
- Explore **druggable targets** (kinases, receptors, transporters)

### For Bioinformaticians
- All intermediate files preserved for custom analysis
- Network data compatible with Cytoscape, Gephi, etc.
- Modular scoring system allows criterion weighting
- Reproducible parameter settings throughout

## 🎯 Clinical Relevance

### Diagnostic Potential
- **74 mitochondrial biomarkers** for patient stratification
- **Complex I/IV signature** for subtype classification
- **Extracellular vesicle proteins** as non-invasive biomarkers

### Therapeutic Opportunities
- **37 druggable targets** including kinases and receptors
- **PARP1 inhibition** potential based on network centrality
- **MCM complex** as DNA replication targets

## 📈 Quality Control

### Technical Validation
- **Excellent data quality**: Sample correlations 0.876-0.977
- **Low missingness**: 7.7% overall, well-handled by MinProb imputation
- **Proper normalization**: Median-centered distributions
- **Robust statistics**: Empirical Bayes moderation for small samples

### Biological Validation
- **Recapitulates known biology**: Mitochondrial complexes affected as expected
- **Literature consistency**: NDUFAF5 identified as top candidate (known Leigh gene)
- **Network validation**: STRING web interface confirms interactions
- **Multi-method convergence**: Enrichment and network analysis tell consistent story

## 🚨 Limitations & Considerations

### Analytical
- Small sample size (n=3 per group) limits statistical power
- Missing value imputation assumes MNAR pattern
- STRING network limited to known interactions

### Biological
- iPSC-derived model may not capture all in vivo complexity
- Focus on proteomics only (integration with transcriptomics recommended)
- Time-point specific (developmental effects may vary)

## 🔮 Future Directions

### Immediate Priorities
1. Validate top candidates via Western blot/qPCR
2. Test mitochondrial function in patient cells
3. Develop biomarker assays for clinical validation

### Long-term Goals
1. Integrate with transcriptomics and metabolomics
2. Develop machine learning classifiers for patient stratification
3. Explore drug repurposing using connectivity mapping

## 📚 Citation

If you use this pipeline, please cite:
```bibtex
@article{inak2021defective,
  title={Defective metabolic programming impairs early neuronal morphogenesis in neural cultures and an organoid model of Leigh syndrome},
  author={Inak, Gizem and Rybak-Wolf, Agnieszka and Lisowski, Pawel and others},
  journal={Nature Communications},
  volume={12},
  number={1},
  pages={1929},
  year={2021},
  publisher={Nature Publishing Group}
}
```

## 🤝 Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🆘 Support

For technical issues:
1. Check the troubleshooting guide in `/docs/troubleshooting.md`
2. Open an issue on GitHub with error logs and system information
3. Contact the maintainers for specific analysis questions

---

**🚀 Ready for experimental validation!** The pipeline has identified high-confidence biomarkers and provides a clear roadmap for wet-lab confirmation and therapeutic development.