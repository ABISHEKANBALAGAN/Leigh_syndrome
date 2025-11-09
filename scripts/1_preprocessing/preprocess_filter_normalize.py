#!/usr/bin/env python3
"""
Proteomics Preprocessing Pipeline:
- Filter low-confidence proteins
- Log2 transformation
- Normalization
- Missing value imputation
- QC figures
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# Create output directory
output_dir = Path('preprocessing_results')
output_dir.mkdir(exist_ok=True)

print("="*70)
print("PROTEOMICS PREPROCESSING PIPELINE")
print("="*70)

# ===================================================================
# 1. LOAD DATA
# ===================================================================

print("\n[Step 1] Loading data...")
lfq_matrix = pd.read_csv('lfq_matrix_with_annotations.csv')
metadata = pd.read_csv('data/metadata.tsv', sep='\t')

# Extract intensity columns (sample names)
sample_cols = [col for col in lfq_matrix.columns 
               if col not in ['Protein_IDs', 'Gene_names', 'Protein_names']]

print(f"✓ Loaded {len(lfq_matrix)} proteins")
print(f"✓ Found {len(sample_cols)} samples: {sample_cols}")

# Separate protein info from intensities
protein_info = lfq_matrix[['Protein_IDs', 'Gene_names', 'Protein_names']].copy()
intensities = lfq_matrix[sample_cols].copy()

print(f"\nIntensity matrix shape: {intensities.shape}")

# ===================================================================
# 2. SIMPLE FILTERING (ALTERNATIVE)
# ===================================================================

print("\n[Step 2] Filtering low-confidence proteins...")

# Convert 0 to NaN
intensities_na = intensities.replace(0, np.nan)

# Store initial count
n_before = len(intensities_na)

# Filter 1: Remove proteins with all missing values
all_missing = intensities_na.isna().all(axis=1)
intensities_filtered = intensities_na[~all_missing].copy()
protein_info_filtered = protein_info[~all_missing].copy()

print(f"✓ Removed {all_missing.sum()} proteins with all missing values")

# Filter 2: Simple filter - require at least 3 valid values total
valid_counts = intensities_filtered.notna().sum(axis=1)
keep_mask = valid_counts >= 3

intensities_filtered = intensities_filtered[keep_mask].copy()
protein_info_filtered = protein_info_filtered[keep_mask].copy()

# Reset indices for both
intensities_filtered = intensities_filtered.reset_index(drop=True)
protein_info_filtered = protein_info_filtered.reset_index(drop=True)

n_after = len(intensities_filtered)
print(f"✓ Filtered to {n_after} proteins (removed {n_before - n_after})")
print(f"  Criterion: ≥3 valid values across all samples")
# ===================================================================
# 3. LOG2 TRANSFORMATION
# ===================================================================

print("\n[Step 3] Log2 transformation...")

# Log2 transform (adding small constant to avoid log(0))
intensities_log2 = np.log2(intensities_filtered + 1)

print(f"✓ Log2 transformation complete")
print(f"  Range: [{intensities_log2.min().min():.2f}, {intensities_log2.max().max():.2f}]")

# Plot: Before vs After transformation
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Before transformation (raw intensities, non-zero only)
for col in sample_cols:
    non_zero = intensities_filtered[col][intensities_filtered[col] > 0]
    axes[0].hist(non_zero, bins=50, alpha=0.3, label=col)
axes[0].set_xlabel('Raw Intensity')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Before Log2 Transformation')
axes[0].set_yscale('log')
axes[0].legend(fontsize=7, ncol=2)

# After transformation
for col in sample_cols:
    non_zero = intensities_log2[col][intensities_log2[col] > 0]
    axes[1].hist(non_zero, bins=50, alpha=0.3, label=col)
axes[1].set_xlabel('Log2 Intensity')
axes[1].set_ylabel('Frequency')
axes[1].set_title('After Log2 Transformation')
axes[1].legend(fontsize=7, ncol=2)

plt.tight_layout()
plt.savefig(output_dir / '01_log2_transformation.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/01_log2_transformation.png")

# ===================================================================
# 4. MISSING VALUE DIAGNOSTICS
# ===================================================================

print("\n[Step 4] Missing value diagnostics...")

# Replace 0 with NaN for missing value analysis
intensities_log2_na = intensities_log2.replace(0, np.nan)

# Calculate missing value statistics
missing_per_protein = intensities_log2_na.isna().sum(axis=1)
missing_per_sample = intensities_log2_na.isna().sum(axis=0)
total_missing_pct = (intensities_log2_na.isna().sum().sum() / intensities_log2_na.size) * 100

print(f"\nMissing value summary:")
print(f"  Total missing: {total_missing_pct:.2f}%")
print(f"\nMissing per sample:")
for sample, count in missing_per_sample.items():
    pct = (count / len(intensities_log2_na)) * 100
    print(f"  {sample}: {count} ({pct:.1f}%)")

# Plot: Missing value patterns
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Missing values per sample
axes[0, 0].bar(range(len(missing_per_sample)), missing_per_sample.values)
axes[0, 0].set_xticks(range(len(missing_per_sample)))
axes[0, 0].set_xticklabels(missing_per_sample.index, rotation=45, ha='right')
axes[0, 0].set_ylabel('Number of Missing Values')
axes[0, 0].set_title('Missing Values per Sample')
axes[0, 0].grid(axis='y', alpha=0.3)

# Plot 2: Missing values per protein histogram
axes[0, 1].hist(missing_per_protein, bins=30, edgecolor='black')
axes[0, 1].set_xlabel('Number of Missing Samples')
axes[0, 1].set_ylabel('Number of Proteins')
axes[0, 1].set_title('Distribution of Missing Values per Protein')
axes[0, 1].grid(axis='y', alpha=0.3)

# Plot 3: Missing value heatmap (sample first 500 proteins for visibility)
n_show = min(500, len(intensities_log2_na))
missing_mask = intensities_log2_na.head(n_show).isna().T
sns.heatmap(missing_mask, cmap='RdYlGn_r', cbar_kws={'label': 'Missing'}, 
            yticklabels=True, xticklabels=False, ax=axes[1, 0])
axes[1, 0].set_title(f'Missing Value Pattern (first {n_show} proteins)')
axes[1, 0].set_xlabel('Protein Index')
axes[1, 0].set_ylabel('Sample')

# Plot 4: Cumulative missing values
sorted_missing = np.sort(missing_per_protein)
cumulative = np.arange(1, len(sorted_missing) + 1) / len(sorted_missing) * 100
axes[1, 1].plot(sorted_missing, cumulative)
axes[1, 1].set_xlabel('Number of Missing Samples')
axes[1, 1].set_ylabel('Cumulative % of Proteins')
axes[1, 1].set_title('Cumulative Distribution of Missing Values')
axes[1, 1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / '02_missing_value_diagnostics.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/02_missing_value_diagnostics.png")

# ===================================================================
# 5. NORMALIZATION
# ===================================================================

print("\n[Step 5] Normalization (Median Centering)...")

# Calculate median of each sample
medians = intensities_log2_na.median()
print("\nSample medians before normalization:")
for sample, med in medians.items():
    print(f"  {sample}: {med:.2f}")

# Median centering normalization
intensities_normalized = intensities_log2_na.sub(medians, axis=1)

print("\n✓ Normalization complete")

# Verify normalization
medians_after = intensities_normalized.median()
print("\nSample medians after normalization (should be ~0):")
for sample, med in medians_after.items():
    print(f"  {sample}: {med:.4f}")

# Plot: Before vs After normalization
fig, axes = plt.subplots(2, 2, figsize=(16, 10))

# Boxplots before normalization
intensities_log2_na.boxplot(ax=axes[0, 0], rot=45)
axes[0, 0].set_title('Before Normalization - Boxplot')
axes[0, 0].set_ylabel('Log2 Intensity')

# Boxplots after normalization
intensities_normalized.boxplot(ax=axes[0, 1], rot=45)
axes[0, 1].set_title('After Normalization - Boxplot')
axes[0, 1].set_ylabel('Normalized Log2 Intensity')

# Density plots before
for col in sample_cols:
    data = intensities_log2_na[col].dropna()
    axes[1, 0].hist(data, bins=50, alpha=0.3, label=col, density=True)
axes[1, 0].set_xlabel('Log2 Intensity')
axes[1, 0].set_ylabel('Density')
axes[1, 0].set_title('Before Normalization - Distribution')
axes[1, 0].legend(fontsize=7, ncol=2)

# Density plots after
for col in sample_cols:
    data = intensities_normalized[col].dropna()
    axes[1, 1].hist(data, bins=50, alpha=0.3, label=col, density=True)
axes[1, 1].set_xlabel('Normalized Log2 Intensity')
axes[1, 1].set_ylabel('Density')
axes[1, 1].set_title('After Normalization - Distribution')
axes[1, 1].legend(fontsize=7, ncol=2)

plt.tight_layout()
plt.savefig(output_dir / '03_normalization_effects.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/03_normalization_effects.png")

# ===================================================================
# 6. IMPUTATION
# ===================================================================

print("\n[Step 6] Missing value imputation...")
print("\nImputation method: MinProb (Left-censored distribution)")

# MinProb imputation: impute from downshifted normal distribution
intensities_imputed = intensities_normalized.copy()

for col in sample_cols:
    col_data = intensities_imputed[col]
    missing_mask = col_data.isna()
    n_missing = missing_mask.sum()
    
    if n_missing > 0:
        # Calculate parameters from observed data
        observed = col_data.dropna()
        mean_obs = observed.mean()
        std_obs = observed.std()
        
        # Impute from left-censored distribution
        # Mean: shifted down by 1.8 standard deviations
        # Std: 30% of original std
        impute_mean = mean_obs - 1.8 * std_obs
        impute_std = 0.3 * std_obs
        
        # Generate random values
        np.random.seed(42)  # For reproducibility
        imputed_values = np.random.normal(impute_mean, impute_std, n_missing)
        
        # Fill missing values
        intensities_imputed.loc[missing_mask, col] = imputed_values
        
        print(f"  {col}: Imputed {n_missing} values")
        print(f"    Mean: {impute_mean:.2f}, Std: {impute_std:.2f}")

print("\n✓ Imputation complete")
print(f"  Remaining NaN: {intensities_imputed.isna().sum().sum()}")

# Plot: Before vs After imputation
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Select first 3 samples for detailed view
samples_to_plot = sample_cols[:3]

for idx, sample in enumerate(samples_to_plot):
    # Before imputation
    axes[0, idx].hist(intensities_normalized[sample].dropna(), 
                     bins=50, alpha=0.7, label='Observed', color='blue')
    axes[0, idx].set_xlabel('Normalized Log2 Intensity')
    axes[0, idx].set_ylabel('Frequency')
    axes[0, idx].set_title(f'{sample} - Before Imputation')
    axes[0, idx].legend()
    
    # After imputation
    # Separate observed and imputed for visualization
    was_missing = intensities_normalized[sample].isna()
    observed_vals = intensities_imputed[sample][~was_missing]
    imputed_vals = intensities_imputed[sample][was_missing]
    
    axes[1, idx].hist(observed_vals, bins=50, alpha=0.7, 
                     label='Observed', color='blue')
    axes[1, idx].hist(imputed_vals, bins=30, alpha=0.7, 
                     label='Imputed', color='red')
    axes[1, idx].set_xlabel('Normalized Log2 Intensity')
    axes[1, idx].set_ylabel('Frequency')
    axes[1, idx].set_title(f'{sample} - After Imputation')
    axes[1, idx].legend()

plt.tight_layout()
plt.savefig(output_dir / '04_imputation_effects.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/04_imputation_effects.png")

# ===================================================================
# 7. SAMPLE CORRELATION & PCA
# ===================================================================

print("\n[Step 7] Sample correlation and PCA...")

# Correlation heatmap
correlation_matrix = intensities_imputed.corr()

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Correlation heatmap
sns.heatmap(correlation_matrix, annot=True, fmt='.2f', cmap='coolwarm',
            center=0, vmin=0, vmax=1, ax=axes[0], 
            cbar_kws={'label': 'Pearson Correlation'})
axes[0].set_title('Sample Correlation (After Preprocessing)')

# PCA
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Transpose (samples as rows)
data_for_pca = intensities_imputed.T

# Standardize
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data_for_pca)

# PCA
pca = PCA(n_components=min(5, len(sample_cols)))
pcs = pca.fit_transform(data_scaled)

# Plot PC1 vs PC2
scatter = axes[1].scatter(pcs[:, 0], pcs[:, 1], s=150, alpha=0.7, c=range(len(sample_cols)), cmap='tab10')
for i, sample in enumerate(sample_cols):
    axes[1].annotate(sample, (pcs[i, 0], pcs[i, 1]), fontsize=9, ha='right')

axes[1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
axes[1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
axes[1].set_title('PCA - Preprocessed Data')
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / '05_correlation_pca.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/05_correlation_pca.png")

# ===================================================================
# 8. SAVE PREPROCESSED DATA
# ===================================================================

print("\n[Step 8] Saving preprocessed data...")

# Combine protein info with preprocessed intensities
preprocessed_full = pd.concat([protein_info_filtered, intensities_imputed], axis=1)

# Save all intermediate and final files
preprocessed_full.to_csv(output_dir / 'preprocessed_data_full.csv', index=False)
print(f"✓ Saved: {output_dir}/preprocessed_data_full.csv")

intensities_imputed.to_csv(output_dir / 'intensities_imputed.csv', index=False)
print(f"✓ Saved: {output_dir}/intensities_imputed.csv")

protein_info_filtered.to_csv(output_dir / 'protein_info_filtered.csv', index=False)
print(f"✓ Saved: {output_dir}/protein_info_filtered.csv")

# Save as pickle for easy loading in Python
import pickle
preprocessing_data = {
    'intensities_imputed': intensities_imputed,
    'intensities_normalized': intensities_normalized,
    'intensities_log2': intensities_log2_na,
    'protein_info': protein_info_filtered,
    'metadata': metadata,
    'sample_cols': sample_cols
}

with open(output_dir / 'preprocessed_data.pkl', 'wb') as f:
    pickle.dump(preprocessing_data, f)
print(f"✓ Saved: {output_dir}/preprocessed_data.pkl")

# ===================================================================
# 9. SUMMARY REPORT
# ===================================================================

print("\n" + "="*70)
print("PREPROCESSING SUMMARY REPORT")
print("="*70)

print(f"\nInput:")
print(f"  Total proteins: {n_before}")
print(f"  Samples: {len(sample_cols)}")

print(f"\nFiltering:")
print(f"  Proteins after filtering: {n_after}")
print(f"  Proteins removed: {n_before - n_after}")
print(f"  Removal rate: {(n_before - n_after) / n_before * 100:.1f}%")

print(f"\nTransformation:")
print(f"  Method: Log2")
print(f"  Range after log2: [{intensities_log2.min().min():.2f}, {intensities_log2.max().max():.2f}]")

print(f"\nNormalization:")
print(f"  Method: Median centering")
print(f"  Median after normalization: ~0 (verified)")

print(f"\nImputation:")
print(f"  Method: MinProb (left-censored distribution)")
print(f"  Missing values before: {intensities_log2_na.isna().sum().sum()}")
print(f"  Missing values after: {intensities_imputed.isna().sum().sum()}")

print(f"\nQuality Control:")
print(f"  Sample correlation range: [{correlation_matrix.values[np.triu_indices_from(correlation_matrix.values, k=1)].min():.3f}, "
      f"{correlation_matrix.values[np.triu_indices_from(correlation_matrix.values, k=1)].max():.3f}]")
print(f"  PC1 explains: {pca.explained_variance_ratio_[0]*100:.1f}% variance")
print(f"  PC2 explains: {pca.explained_variance_ratio_[1]*100:.1f}% variance")

print(f"\nOutput Files:")
print(f"  - preprocessed_data_full.csv        (Complete data)")
print(f"  - intensities_imputed.csv           (Matrix only)")
print(f"  - protein_info_filtered.csv         (Annotations)")
print(f"  - preprocessed_data.pkl             (Python pickle)")

print(f"\nQC Figures:")
print(f"  - 01_log2_transformation.png")
print(f"  - 02_missing_value_diagnostics.png")
print(f"  - 03_normalization_effects.png")
print(f"  - 04_imputation_effects.png")
print(f"  - 05_correlation_pca.png")

print("\n" + "="*70)
print("PREPROCESSING COMPLETE!")
print("="*70)
print(f"\nNext step: Differential abundance analysis")
print(f"  Use file: {output_dir}/preprocessed_data.pkl")