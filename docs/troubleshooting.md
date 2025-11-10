I'll integrate the Docker troubleshooting content into your existing troubleshooting guide. Here's the updated comprehensive version:

# Troubleshooting Guide

## Common Issues and Solutions

### 🐳 Docker Issues

#### Docker Build Failures
```bash
# Problem: Docker build fails
# Solution:
# 1. Ensure Docker is running
docker --version

# 2. Check disk space
docker system df

# 3. Try rebuilding with cache disabled
docker build --no-cache -t leigh_pipeline .

# 4. Check Docker daemon status
docker info
```

#### Container File Access Issues
```bash
# Problem: Container can't access files
# Solution:
# 1. Check file permissions in container
docker run -it leigh_pipeline ls -la /app/

# 2. Ensure data files are copied correctly
docker run -it leigh_pipeline ls /app/data/

# 3. Run with volume mount for direct file access
docker run -v $(pwd):/app -it leigh_pipeline

# 4. Verify file ownership
docker run -it leigh_pipeline chown -R 1000:1000 /app/data
```

#### Docker Memory Issues
```bash
# For large datasets, increase Docker resources
docker run -it --memory=8g --cpus=4 leigh_pipeline

# Check current container resource usage
docker stats

# Clean up unused containers and images
docker system prune
```

### 🐍 Python Script Issues

#### ModuleNotFoundError
```bash
# If you get "No module named 'pandas'"
pip install pandas numpy matplotlib seaborn scipy statsmodels scikit-learn requests

# Or using conda:
conda install pandas numpy matplotlib seaborn scipy statsmodels scikit-learn requests

# In Docker container, check installed packages:
docker run -it leigh_pipeline pip list

# Check Python path:
docker run -it leigh_pipeline python -c "import sys; print(sys.path)"
```

#### FileNotFoundError
```bash
# Ensure you're in the correct directory
pwd  # Leigh_syndrome

# Check if data files exist in container
docker run -it leigh_pipeline ls -la /app/data/
# Should show: proteinGroups_clean.txt, metadata.tsv

# Verify file encoding:
docker run -it leigh_pipeline file -i /app/data/proteinGroups_clean.txt
```

#### Memory Errors
```bash
# For large datasets, increase memory
python -W ignore scripts/1_preprocessing/preprocess_filter_normalize.py

# Or process in chunks by modifying the script:
# Add: pandas.read_csv(..., chunksize=1000)

# In Docker, increase container memory:
docker run -it --memory=8g leigh_pipeline

# Clear memory explicitly in scripts:
import gc
gc.collect()
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

# For Docker environments:
BiocManager::install("STRINGdb", configure.args="--disable-threading")

# Update Bioconductor if needed:
BiocManager::install(version = "3.17")
```

#### STRING Mapping Failures
```R
# If gene mapping fails, check the input format:
# Ensure genes_df is a data.frame (not tibble)
genes_df <- as.data.frame(genes_df)
genes_df$Gene_names <- as.character(genes_df$Gene_names)

# Check for non-standard gene symbols
unique(nchar(genes_df$Gene_names))  # Should be reasonable lengths

# Handle network timeouts:
options(timeout = 400)  # Increase timeout
```

#### Network Timeout Issues
```R
# Problem: Network timeout in STRINGdb
# Solution:
# 1. Increase timeout significantly
options(timeout = 600)

# 2. Check internet connectivity in container
system("ping -c 3 biit.cs.ut.ee")

# 3. Retry with smaller gene sets
MAX_GENES <- 100  # Reduce from 200 if needed

# 4. Lower confidence threshold for faster results
SCORE_THRESHOLD <- 300
```

### 🔧 Data Preparation Issues

#### proteinGroups_clean.txt Creation
```bash
# If the grep command fails, check the original file structure:
head -n 5 data/proteinGroups.txt
# Look for "CON__" and "REV__" patterns

# Alternative filtering approach:
awk '! /CON__|REV__/ {print}' data/proteinGroups.txt > data/proteinGroups_clean.txt

# In Docker, ensure files are accessible:
docker run -v $(pwd)/data:/app/data -it leigh_pipeline head -n 5 /app/data/proteinGroups.txt
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

# Check data format and column names:
print([col for col in df.columns if 'LFQ' in col])
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

# 4. Check mapping rate and adjust gene set
mapping_rate <- sum(!is.na(mapped$STRING_id)) / nrow(mapped)
print(paste("Mapping rate:", round(mapping_rate * 100, 1), "%"))
```

#### Cytoscape File Issues
**Problem**: Cytoscape files not loading properly
```R
# Solutions:
# 1. Check file format (should be tab-separated)
write.table(nodes, "nodes.txt", sep="\t", quote=FALSE, row.names=FALSE)

# 2. Verify column headers match expected format
# 3. Import network first, then attributes in Cytoscape
# 4. Check for special characters in gene names
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

# In Docker on Windows, use PowerShell syntax:
docker run -v ${PWD}:/app -it leigh_pipeline
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

# Optimize Docker performance:
docker run -it --memory=8g --cpus=4 leigh_pipeline python scripts/your_script.py
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

#### Slow Network Analysis Optimization
```R
# Reduce gene set size for faster processing
MAX_GENES <- 100  # Instead of 200 for initial testing

# Lower confidence threshold for more results
SCORE_THRESHOLD <- 300  # Instead of 400

# Process in stages for large analyses
batch_size <- 50
for(i in seq(1, nrow(genes_df), batch_size)) {
    batch <- genes_df[i:min(i+batch_size-1, nrow(genes_df)), ]
    # Process batch
}
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

# In Docker, capture all output:
docker run -it leigh_pipeline python script.py 2>&1 | tee log.txt
```

#### Check Intermediate Files
```bash
# Verify each step produces expected outputs
ls -la preprocessing_results/  # Should have 9 files after Day 1
ls -la differential_results/   # Should have 9 files after Day 2

# In Docker, check container file system:
docker run -it leigh_pipeline ls -la /app/preprocessing_results/
```

#### Validate Data Quality
```python
# Add data quality checks to scripts
def validate_data(df, stage):
    assert not df.isnull().all().any(), f"All-NA columns in {stage}"
    assert df.shape[0] > 0, f"Empty dataframe in {stage}"
    print(f"Data validation passed for {stage}: {df.shape}")
    return True

# Check for NaN/Inf values explicitly:
import numpy as np
print(f"NaN values: {df.isna().sum().sum()}")
print(f"Inf values: {np.isinf(df.select_dtypes(include=[np.number])).sum().sum()}")
```

### 📋 Common Error Messages and Solutions

#### "KeyError: 'LFQ intensity'"
**Cause**: Column names don't match expected pattern
**Solution**: Check MaxQuant output format
```python
# List all columns to find correct LFQ columns
print([col for col in pg.columns if 'LFQ' in col])

# Alternative column patterns
lfq_cols = [col for col in pg.columns if 'intensity' in col.lower()]
```

#### "HTTP Error 429: Too Many Requests"
**Cause**: g:Profiler API rate limiting
**Solution**: Add delays between requests
```python
import time
time.sleep(2)  # Wait 2 seconds between requests

# Or use exponential backoff
import random
time.sleep(random.uniform(1, 3))
```

#### "Network is too large for STRING"
**Cause**: Too many input genes (>200)
**Solution**: Filter to top significant genes
```R
# Keep only top genes by significance
top_genes <- head(genes_df[order(genes_df$adj.P.Val), ], 150)

# Or use more stringent significance cutoff
significant_genes <- genes_df[genes_df$adj.P.Val < 0.01, ]
```

#### "Error in string_db$map()"
**Cause**: Gene symbol format issues or STRINGdb version mismatch
**Solution**: 
```R
# Check gene symbols are official HGNC symbols
# Verify STRINGdb version matches current database
# Try with a small test set first
test_genes <- head(genes_df, 10)
mapped_test <- string_db$map(test_genes, "Gene_names")
```

### 🆘 Getting Help

If you encounter issues not covered here:

1. **Check the logs**: All scripts produce detailed output
2. **Verify file permissions**: `ls -la` to check read/write permissions
3. **Check system resources**: `free -h` and `df -h` for memory/disk space
4. **Create a minimal example**: Reproduce with small test dataset
5. **Open an issue**: Include error message and system information

#### Docker-Specific Help Steps
```bash
# 1. Check Docker container status
docker ps -a

# 2. Examine container logs
docker logs <container_id>

# 3. Test with minimal data subset
docker run -it leigh_pipeline python -c "import pandas as pd; print('Dependencies OK')"

# 4. Verify data files in container
docker run -v $(pwd)/data:/app/data -it leigh_pipeline ls -la /app/data/

# 5. Check resource allocation
docker system info
```

### Example Debugging Session
```bash
# Start from scratch with Docker
cd Leigh_syndrome

# Build fresh image
docker build --no-cache -t leigh_pipeline .

# Verify data files are accessible
docker run -v $(pwd)/data:/app/data -it leigh_pipeline ls -la /app/data/
# Should see: proteinGroups_clean.txt, metadata.tsv

# Run step by step with verbose output
docker run -v $(pwd)/data:/app/data -it leigh_pipeline python -v /app/scripts/0_processing/preprocess.py

# Check intermediate results
docker run -v $(pwd):/app -it leigh_pipeline head /app/preprocessing_results/preprocessed_data_full.csv

# Continue to next step with increased resources
docker run -v $(pwd):/app --memory=8g --cpus=4 -it leigh_pipeline python /app/scripts/1_preprocessing/preprocess_filter_normalize.py
```

## Process in Stages for Better Debugging
```bash
# Stage 1: Data preparation only
docker run -v $(pwd):/app -it leigh_pipeline python /app/scripts/0_processing/preprocess.py

# Stage 2: Preprocessing
docker run -v $(pwd):/app -it leigh_pipeline python /app/scripts/1_preprocessing/preprocessing_pipeline.py

# Stage 3: Differential analysis
docker run -v $(pwd):/app -it leigh_pipeline python /app/scripts/2_differential_analysis/differential_analysis.py

# Stage 4: Network analysis (can be memory intensive)
docker run -v $(pwd):/app --memory=8g -it leigh_pipeline Rscript /app/scripts/3_Functional_enrichment_and_Biological_annotation/string_ppi.R
```

This comprehensive troubleshooting guide now covers both standard execution and Docker-specific issues. The key is to run each step sequentially, verify outputs, and use the Docker volume mounts for easy file access and debugging.