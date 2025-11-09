#!/usr/bin/env python3
"""
Import MaxQuant proteinGroups_clean.txt and extract LFQ intensity matrix
"""

import pandas as pd
import numpy as np

# ===================================================================
# 1. LOAD DATA
# ===================================================================

print("Loading proteinGroups_clean.txt...")
pg = pd.read_csv('data/proteinGroups_clean.txt', sep='\t', low_memory=False)
print(f"✓ Loaded {len(pg)} proteins x {len(pg.columns)} columns")

print("\nLoading metadata.tsv...")
metadata = pd.read_csv('data/metadata.tsv', sep='\t')
print(f"✓ Loaded {len(metadata)} samples")
print("\nMetadata:")
print(metadata)

# ===================================================================
# 2. EXTRACT LFQ INTENSITY COLUMNS
# ===================================================================

# Get all LFQ intensity columns
lfq_cols = [col for col in pg.columns if col.startswith('LFQ intensity ')]
print(f"\n✓ Found {len(lfq_cols)} LFQ intensity columns:")
print(lfq_cols)

# Extract LFQ matrix
lfq_matrix = pg[lfq_cols].copy()

# Get sample names (remove 'LFQ intensity ' prefix)
sample_names = [col.replace('LFQ intensity ', '') for col in lfq_cols]
lfq_matrix.columns = sample_names

print(f"\n✓ Extracted LFQ matrix: {lfq_matrix.shape}")
print(f"   Samples: {sample_names}")

# ===================================================================
# 3. ADD PROTEIN ANNOTATIONS
# ===================================================================

# Add protein identifiers and gene names as index/columns
lfq_matrix.insert(0, 'Protein_IDs', pg['Protein IDs'])
lfq_matrix.insert(1, 'Gene_names', pg['Gene names'])
lfq_matrix.insert(2, 'Protein_names', pg['Protein names'])

print(f"\n✓ Added protein annotations")

# ===================================================================
# 4. DATA SUMMARY
# ===================================================================

print("\n" + "="*60)
print("DATA SUMMARY")
print("="*60)

# Matrix dimensions
print(f"Proteins: {len(lfq_matrix)}")
print(f"Samples: {len(sample_names)}")

# Missing values (0 or NA in LFQ data)
lfq_data_only = lfq_matrix[sample_names]
lfq_data_only_numeric = lfq_data_only.replace(0, np.nan)
missing_per_sample = lfq_data_only_numeric.isna().sum()
missing_pct_per_sample = (missing_per_sample / len(lfq_matrix)) * 100

print(f"\nMissing values per sample:")
for sample, count, pct in zip(sample_names, missing_per_sample, missing_pct_per_sample):
    print(f"  {sample}: {count} ({pct:.1f}%)")

# Proteins with all zeros
all_zero = (lfq_data_only == 0).all(axis=1).sum()
print(f"\nProteins with all zeros: {all_zero}")

# Value ranges
print(f"\nIntensity ranges (non-zero):")
for sample in sample_names:
    non_zero = lfq_matrix[sample][lfq_matrix[sample] > 0]
    if len(non_zero) > 0:
        print(f"  {sample}: {non_zero.min():.2e} - {non_zero.max():.2e}")

# ===================================================================
# 5. SAVE OUTPUTS
# ===================================================================

print("\n" + "="*60)
print("SAVING FILES")
print("="*60)

# Save complete LFQ matrix with annotations
lfq_matrix.to_csv('lfq_matrix_with_annotations.csv', index=False)
print("✓ Saved: lfq_matrix_with_annotations.csv")

# Save intensity-only matrix (for analysis)
lfq_data_only.to_csv('lfq_matrix_intensities_only.csv', index=False)
print("✓ Saved: lfq_matrix_intensities_only.csv")

# Save protein info separately
protein_info = pg[['Protein IDs', 'Majority protein IDs', 'Gene names', 
                   'Protein names', 'Peptides', 'Sequence coverage [%]',
                   'Mol. weight [kDa]']].copy()
protein_info.to_csv('protein_info.csv', index=False)
print("✓ Saved: protein_info.csv")

# ===================================================================
# 6. PREVIEW
# ===================================================================

print("\n" + "="*60)
print("PREVIEW: First 5 rows of LFQ matrix")
print("="*60)
print(lfq_matrix.head())

print("\n" + "="*60)
print("COMPLETE!")
print("="*60)
print("\nOutput files:")
print("  1. lfq_matrix_with_annotations.csv  - Full matrix with protein info")
print("  2. lfq_matrix_intensities_only.csv  - Intensity values only")
print("  3. protein_info.csv                 - Protein annotations")
print("\nNext step: Use these files for preprocessing and analysis")