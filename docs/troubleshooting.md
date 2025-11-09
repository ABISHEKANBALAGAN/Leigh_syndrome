I'll create a comprehensive troubleshooting guide for you. Here's the content for `/docs/troubleshooting.md`:

# Troubleshooting Guide

## Common Issues and Solutions

### 🐍 Python Script Issues

#### ModuleNotFoundError
```bash
# If you get "No module named 'pandas'"
pip install pandas numpy matplotlib seaborn scipy statsmodels scikit-learn requests

# Or using conda:
conda install pandas numpy matplotlib seaborn scipy statsmodels scikit-learn requests
```

#### FileNotFoundError
```bash
# Ensure you're in the correct directory
pwd  # Leigh_syndrome

# Check if data files exist
ls -la data/
# Should show: proteinGroups_clean.txt, metadata.tsv
```

#### Memory Errors
```bash
# For large datasets, increase memory
python -W ignore scripts/1_preprocessing/preprocess_filter_normalize.py

# Or process in chunks by modifying the script:
# Add: pandas.read_csv(..., chunksize=1000)
```

### R Script Issues

#### STRINGdb Installation Problems
```R
# If BiocManager fails:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install("STRINGdb", version = "3.16")

# Alternative installation:
install.packages("remotes")
remotes::install_bioc("STRINGdb")
```

#### STRING Mapping Failures
```R
# If gene mapping fails, check the input format:
# Ensure genes_df is a data.frame (not tibble)
genes_df <- as.data.frame(genes_df)
genes_df$Gene_names <- as.character(genes_df$Gene_names)

# Check for non-standard gene symbols
unique(nchar(genes_df$Gene_names))  # Should be reasonable lengths
```

### 🔧 Data Preparation Issues

#### proteinGroups_clean.txt Creation
```bash
# If the grep command fails, check the original file structure:
head -n 5 data/proteinGroups.txt
# Look for "CON__" and "REV__" patterns

# Alternative filtering approach:
awk '! /CON__|REV__/ {print}' data/proteinGroups.txt > data/proteinGroups_clean.txt
```

#### Metadata File Creation
```bash
# If awk script fails, check summary.txt format:
head -n 10 data/summary.txt

# Manual metadata creation example:
echo -e "RawFile\tExperiment\tFraction\tGroup\tCondition\tReplicate"
echo -e "3P_1_D1\t3P_1\t1\t3P\tPatient\t1"
echo -e "SIC_1_D1\tSIC_1\t1\tSIC\tControl\t1"
```

### 📊 Analysis-Specific Issues

#### Preprocessing Failures
**Problem**: Log2 transformation errors with negative values
```python
# Solution: Ensure no negative values before log2
intensities_filtered = intensities_filtered.clip(lower=0)
intensities_log2 = np.log2(intensities_filtered + 1)
```

**Problem**: Normalization issues with many missing values
```python
# Solution: Use more stringent filtering
valid_counts = intensities_filtered.notna().sum(axis=1)
keep_mask = valid_counts >= 4  # Require more valid values
```

#### Differential Analysis Warnings
**Problem**: "Division by zero" in moderated t-test
```python
# Solution: Add small constant to variances
pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / df_residual
pooled_var = np.where(pooled_var == 0, 1e-10, pooled_var)
```

#### Network Analysis Problems
**Problem**: No interactions found in STRING
```R
# Solutions:
# 1. Lower the confidence threshold
string_db <- STRINGdb$new(score_threshold = 300)  # Lower from 400

# 2. Check gene symbol mapping
mapped <- string_db$map(genes_df, "Gene_names", removeUnmappedRows = FALSE)
print(mapped[is.na(mapped$STRING_id), ])  # Show unmapped genes

# 3. Use different identifier types
# Try ENSEMBL IDs or UniProt accessions if available
```

### 🖥️ System-Specific Issues

#### Linux/Mac Path Issues
```bash
# If scripts aren't executable:
chmod +x scripts/*/*.py
chmod +x scripts/*/*.R

# Python shebang issues:
which python3
# Should show: /usr/bin/python3 or /home/username/anaconda3/bin/python3
```

#### Windows Compatibility
```bash
# For Windows users, use Python from Anaconda prompt
# Replace forward slashes with backslashes in paths
python scripts\1_preprocessing\preprocess_filter_normalize.py

# For R on Windows, use Rscript from command line
Rscript scripts\3_Functional_enrichment_and_Biological_annotation\string_ppi.R
```

### 📈 Performance Optimization

#### Large Dataset Handling
```python
# Process in chunks for memory efficiency
chunk_size = 1000
for chunk in pd.read_csv('large_file.csv', chunksize=chunk_size):
    process_chunk(chunk)

# Use data types to reduce memory
dtypes = {'Protein_IDs': 'string', 'Gene_names': 'string'}
lfq_matrix = pd.read_csv('file.csv', dtype=dtypes)
```

#### Parallel Processing
```python
# For computationally intensive steps
from multiprocessing import Pool

def process_protein(protein_data):
    return some_computation(protein_data)

with Pool(processes=4) as pool:
    results = pool.map(process_protein, protein_list)
```

### 🔍 Debugging Tips

#### Enable Verbose Output
```python
# Add debug prints to scripts
import logging
logging.basicConfig(level=logging.DEBUG)

# Or use print statements at key steps
print(f"Processing {len(proteins)} proteins")
print(f"Sample medians: {medians}")
```

#### Check Intermediate Files
```bash
# Verify each step produces expected outputs
ls -la preprocessing_results/  # Should have 9 files after Day 1
ls -la differential_results/   # Should have 9 files after Day 2
```

#### Validate Data Quality
```python
# Add data quality checks to scripts
def validate_data(df, stage):
    assert not df.isnull().all().any(), f"All-NA columns in {stage}"
    assert df.shape[0] > 0, f"Empty dataframe in {stage}"
    return True
```

### 📋 Common Error Messages and Solutions

#### "KeyError: 'LFQ intensity'"
**Cause**: Column names don't match expected pattern
**Solution**: Check MaxQuant output format
```python
# List all columns to find correct LFQ columns
print([col for col in pg.columns if 'LFQ' in col])
```

#### "HTTP Error 429: Too Many Requests"
**Cause**: g:Profiler API rate limiting
**Solution**: Add delays between requests
```python
import time
time.sleep(2)  # Wait 2 seconds between requests
```

#### "Network is too large for STRING"
**Cause**: Too many input genes (>200)
**Solution**: Filter to top significant genes
```R
# Keep only top genes by significance
top_genes <- head(genes_df[order(genes_df$adj.P.Val), ], 150)
```

### 🆘 Getting Help

If you encounter issues not covered here:

1. **Check the logs**: All scripts produce detailed output
2. **Verify file permissions**: `ls -la` to check read/write permissions
3. **Check system resources**: `free -h` and `df -h` for memory/disk space
4. **Create a minimal example**: Reproduce with small test dataset
5. **Open an issue**: Include error message and system information

### Example Debugging Session
```bash
# Start from scratch
cd Leigh_syndrome

# Verify data files
ls -la data/
# Should see: proteinGroups_clean.txt, metadata.tsv

# Run step by step with verbose output
python -v scripts/0_processing/preprocess.py

# Check intermediate results
head preprocessing_results/preprocessed_data_full.csv

# Continue to next step
python -v scripts/1_preprocessing/preprocess_filter_normalize.py
```

This troubleshooting guide should help resolve most common issues with the pipeline. The key is to run each step sequentially and verify the outputs before proceeding to the next step.
