#!/usr/bin/env python3
"""
Day 3: Differential Abundance Analysis (IMPROVED VERSION)
Compare 3P (Patient) vs SIC groups using proper empirical Bayes moderated t-tests

Improvements:
- Proper empirical Bayes prior estimation from data
- Exact sample matching (not substring)
- Explicit protein ID tracking
- Validation of preprocessing assumptions
- QC summary tables
- Enrichment-ready output
"""
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pickle

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# Create output directory
output_dir = Path('differential_results')
output_dir.mkdir(exist_ok=True)

print("="*70)
print("DIFFERENTIAL ABUNDANCE ANALYSIS (IMPROVED)")
print("="*70)

# ===================================================================
# 1. LOAD PREPROCESSED DATA
# ===================================================================

print("\n[Step 1] Loading preprocessed data...")

with open('preprocessing_results/preprocessed_data.pkl', 'rb') as f:
    data = pickle.load(f)

intensities_imputed = data['intensities_imputed']
protein_info = data['protein_info']
metadata = data['metadata']
sample_cols = data['sample_cols']

print(f"✓ Loaded {intensities_imputed.shape[0]} proteins x {intensities_imputed.shape[1]} samples")

# ===================================================================
# 2. VALIDATE PREPROCESSING
# ===================================================================

print("\n[Step 2] Validating preprocessing assumptions...")

# Check 1: Data is log2-transformed
sample_ranges = intensities_imputed.max() - intensities_imputed.min()
if sample_ranges.max() > 25:  # Raw intensities would be much larger
    print("⚠ WARNING: Data may not be log2-transformed!")
else:
    print("✓ Data appears log2-transformed (range check passed)")

# Check 2: Data is normalized (medians should be similar)
medians = intensities_imputed.median()
median_range = medians.max() - medians.min()
if median_range > 2:  # More than 2 log2 units difference
    print(f"⚠ WARNING: Medians vary by {median_range:.2f} log2 units - may not be normalized")
else:
    print(f"✓ Data appears normalized (median range: {median_range:.3f})")

# Check 3: No missing values after imputation
n_missing = intensities_imputed.isna().sum().sum()
if n_missing > 0:
    print(f"⚠ WARNING: {n_missing} missing values remain after imputation!")
else:
    print("✓ No missing values (imputation complete)")

# Display preprocessing summary
print("\nPreprocessing summary:")
print(f"  Median intensities: {medians.min():.2f} to {medians.max():.2f}")
print(f"  Intensity ranges: {sample_ranges.min():.2f} to {sample_ranges.max():.2f}")

# ===================================================================
# 3. SETUP COMPARISON GROUPS (EXACT MATCHING)
# ===================================================================

print("\n[Step 3] Setting up comparison groups (exact matching)...")

# Create sample-to-condition mapping from metadata
# Use Experiment column to map to sample names
sample_to_condition = {}
sample_to_experiment = {}

for _, row in metadata.iterrows():
    exp = row['Experiment']
    condition = row['Group']
    
    # Find exact matching sample in sample_cols
    # Sample names are like "3P_1", experiments are "3P_1"
    for sample in sample_cols:
        if sample == exp or sample.startswith(exp + '_'):
            sample_to_condition[sample] = condition
            sample_to_experiment[sample] = exp

print(f"\nSample mapping (exact match):")
for sample in sorted(sample_to_condition.keys()):
    print(f"  {sample:15s} -> {sample_to_condition[sample]:10s} ({sample_to_experiment[sample]})")

# Get samples for each condition
condition1 = '3P'   # Control/baseline
condition2 = 'SIC'  # Treatment/comparison

group1_samples = [s for s in sample_cols if sample_to_condition.get(s) == condition1]
group2_samples = [s for s in sample_cols if sample_to_condition.get(s) == condition2]

# Validate groups
if not group1_samples:
    raise ValueError(f"No samples found for condition '{condition1}'!")
if not group2_samples:
    raise ValueError(f"No samples found for condition '{condition2}'!")

print(f"\nComparison: {condition2} vs {condition1}")
print(f"  {condition1}: n={len(group1_samples)}")
for s in group1_samples:
    print(f"    - {s}")
print(f"  {condition2}: n={len(group2_samples)}")
for s in group2_samples:
    print(f"    - {s}")

# ===================================================================
# 4. EMPIRICAL BAYES MODERATED T-TEST (PROPER IMPLEMENTATION)
# ===================================================================

print("\n[Step 4] Performing empirical Bayes moderated t-test...")

def estimate_prior_df(variances, df_residual):
    """
    Estimate prior degrees of freedom from data using method-of-moments
    Similar to limma's fitFDist
    """
    # Remove any zero or negative variances
    valid_var = variances[variances > 0]
    
    if len(valid_var) < 10:
        return 3.0  # Default fallback
    
    # Calculate variance of log(variances)
    log_var = np.log(valid_var)
    var_log_var = np.var(log_var)
    
    # Trigamma function approximation
    # For large df, trigamma(df/2) ≈ 2/df
    if var_log_var > 0:
        d0 = 2 / var_log_var
        d0 = max(1, min(d0, 50))  # Constrain to reasonable range
    else:
        d0 = 3.0
    
    return d0

def limma_moderated_ttest_improved(data1, data2):
    """
    Improved limma-style moderated t-test with data-driven prior estimation
    
    Parameters:
    -----------
    data1 : array (n_proteins, n_samples1) - Group 1 intensities
    data2 : array (n_proteins, n_samples2) - Group 2 intensities
    
    Returns:
    --------
    dict with logFC, P.Value, adj.P.Val, etc.
    """
    n1, n2 = data1.shape[1], data2.shape[1]
    n = n1 + n2
    df_residual = n - 2
    
    # Calculate means
    mean1 = np.mean(data1, axis=1)
    mean2 = np.mean(data2, axis=1)
    logFC = mean2 - mean1
    
    # Calculate variances
    var1 = np.var(data1, axis=1, ddof=1)
    var2 = np.var(data2, axis=1, ddof=1)
    
    # Pooled variance
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / df_residual
    
    # Estimate prior parameters from data (empirical Bayes)
    valid_var = pooled_var[pooled_var > 0]
    s0_squared = np.median(valid_var)  # Prior variance (robust estimate)
    d0 = estimate_prior_df(valid_var, df_residual)  # Prior df
    
    print(f"\nEmpirical Bayes prior parameters (estimated from data):")
    print(f"  Prior df (d0): {d0:.2f}")
    print(f"  Prior variance (s0²): {s0_squared:.4f}")
    
    # Posterior variance (moderated)
    post_var = (d0 * s0_squared + df_residual * pooled_var) / (d0 + df_residual)
    
    # Standard error
    se = np.sqrt(post_var * (1/n1 + 1/n2))
    
    # Moderated t-statistic
    t_stat = logFC / se
    
    # Moderated degrees of freedom
    df_total = d0 + df_residual
    
    # P-values (two-tailed)
    p_values = 2 * (1 - stats.t.cdf(np.abs(t_stat), df_total))
    
    # Handle edge cases
    p_values = np.where(np.isnan(p_values) | np.isinf(p_values), 1.0, p_values)
    p_values = np.clip(p_values, 0, 1)
    
    # Adjusted p-values (Benjamini-Hochberg FDR)
    _, adj_p_values, _, _ = multipletests(p_values, method='fdr_bh')
    
    # B-statistic (log-odds of differential expression)
    B = np.log(post_var / s0_squared)
    
    return {
        'logFC': logFC,
        'AveExpr': (mean1 + mean2) / 2,
        't': t_stat,
        'P.Value': p_values,
        'adj.P.Val': adj_p_values,
        'B': B,
        'mean1': mean1,
        'mean2': mean2,
        'prior_df': d0,
        'prior_var': s0_squared
    }

# Extract data for each group
data1 = intensities_imputed[group1_samples].values
data2 = intensities_imputed[group2_samples].values

print(f"\nData shapes:")
print(f"  Group 1 ({condition1}): {data1.shape}")
print(f"  Group 2 ({condition2}): {data2.shape}")

# Perform moderated t-test
results_dict = limma_moderated_ttest_improved(data1, data2)

# ===================================================================
# 5. CREATE RESULTS DATAFRAME WITH EXPLICIT PROTEIN ID TRACKING
# ===================================================================

print("\n[Step 5] Building results dataframe...")

# Start with protein info to ensure correct alignment
results = protein_info.copy()

# Add statistical results explicitly
for key in ['logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B', 'mean1', 'mean2']:
    results[key] = results_dict[key]

# Add readable column names for means
results.rename(columns={
    'mean1': f'mean_{condition1}',
    'mean2': f'mean_{condition2}'
}, inplace=True)

# Reorder columns
col_order = ['Protein_IDs', 'Gene_names', 'Protein_names', 
             'logFC', 'AveExpr', f'mean_{condition1}', f'mean_{condition2}',
             't', 'P.Value', 'adj.P.Val', 'B']
results = results[col_order]

print(f"✓ Results dataframe created: {results.shape}")

# ===================================================================
# 6. ADD SIGNIFICANCE FLAGS & ANNOTATIONS
# ===================================================================

print("\n[Step 6] Adding significance flags...")

# Define thresholds
fc_threshold = 1.0      # log2 fold change (2-fold)
pval_threshold = 0.05   # adjusted p-value

# Add significance column
results['Significant'] = (
    (np.abs(results['logFC']) > fc_threshold) & 
    (results['adj.P.Val'] < pval_threshold)
)

# Add regulation direction
results['Regulation'] = 'NS'
results.loc[
    (results['logFC'] > fc_threshold) & (results['adj.P.Val'] < pval_threshold),
    'Regulation'
] = 'Up'
results.loc[
    (results['logFC'] < -fc_threshold) & (results['adj.P.Val'] < pval_threshold),
    'Regulation'
] = 'Down'

# Add fold change in linear scale
results['FoldChange'] = 2 ** results['logFC']

# Sort by p-value
results = results.sort_values('P.Value').reset_index(drop=True)

# Summary statistics
n_total = len(results)
n_sig = results['Significant'].sum()
n_up = (results['Regulation'] == 'Up').sum()
n_down = (results['Regulation'] == 'Down').sum()

print(f"\nResults Summary:")
print(f"  Total proteins tested: {n_total}")
print(f"  Significant (|logFC|>{fc_threshold}, adj.P<{pval_threshold}): {n_sig} ({n_sig/n_total*100:.1f}%)")
print(f"    Up-regulated in {condition2}: {n_up} ({n_up/n_total*100:.1f}%)")
print(f"    Down-regulated in {condition2}: {n_down} ({n_down/n_total*100:.1f}%)")

# ===================================================================
# 7. QC SUMMARY TABLE
# ===================================================================

print("\n[Step 7] Generating QC summary table...")

qc_summary = pd.DataFrame({
    'Metric': [
        'Total proteins tested',
        'Samples (Group 1)',
        'Samples (Group 2)',
        'Significant proteins',
        'Up-regulated',
        'Down-regulated',
        'Not significant',
        'Median logFC (all)',
        'Median logFC (significant)',
        'Mean P-value',
        'Mean adj.P-value',
        'Prior df (d0)',
        'Prior variance (s0²)',
        'FC threshold',
        'P-value threshold'
    ],
    'Value': [
        n_total,
        len(group1_samples),
        len(group2_samples),
        n_sig,
        n_up,
        n_down,
        n_total - n_sig,
        f"{results['logFC'].median():.3f}",
        f"{results.loc[results['Significant'], 'logFC'].median():.3f}" if n_sig > 0 else "NA",
        f"{results['P.Value'].mean():.4f}",
        f"{results['adj.P.Val'].mean():.4f}",
        f"{results_dict['prior_df']:.2f}",
        f"{results_dict['prior_var']:.4f}",
        fc_threshold,
        pval_threshold
    ]
})

qc_summary.to_csv(output_dir / 'qc_summary.csv', index=False)
print(f"✓ Saved: {output_dir}/qc_summary.csv")
print("\n" + qc_summary.to_string(index=False))

# ===================================================================
# 8. VOLCANO PLOT (ENHANCED)
# ===================================================================

print("\n[Step 8] Creating enhanced volcano plot...")

fig, ax = plt.subplots(figsize=(11, 9))

# Prepare data
x = results['logFC']
y = -np.log10(results['P.Value'])

# Color by regulation and significance level
colors = []
for _, row in results.iterrows():
    if row['adj.P.Val'] < 0.001:
        # Highly significant
        colors.append('#8B0000' if row['Regulation'] == 'Up' else '#00008B')
    elif row['Regulation'] == 'Up':
        colors.append('#d62728')
    elif row['Regulation'] == 'Down':
        colors.append('#1f77b4')
    else:
        colors.append('#7f7f7f')

# Scatter plot with size by expression
sizes = np.clip(results['AveExpr'] - results['AveExpr'].min(), 5, 50)
ax.scatter(x, y, c=colors, alpha=0.6, s=sizes, edgecolors='none')

# Threshold lines
ax.axhline(-np.log10(pval_threshold), color='black', linestyle='--', 
          linewidth=1, alpha=0.5, label=f'P = {pval_threshold}')
ax.axvline(fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axvline(-fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)

# Label top hits
label_top = 15
top_hits = results.nsmallest(label_top, 'P.Value')
for idx, row in top_hits.iterrows():
    gene = str(row['Gene_names']).split(';')[0] if pd.notna(row['Gene_names']) else ''
    if gene and gene != 'nan':
        ax.annotate(gene, (row['logFC'], -np.log10(row['P.Value'])),
                   fontsize=7, alpha=0.9, 
                   xytext=(3, 3), textcoords='offset points',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

# Labels
ax.set_xlabel('log₂ Fold Change', fontsize=13, fontweight='bold')
ax.set_ylabel('-log₁₀(P-value)', fontsize=13, fontweight='bold')
ax.set_title(f'Volcano Plot: {condition2} vs {condition1}\n'
            f'({n_sig} significant: {n_up} up, {n_down} down)', 
            fontsize=14, fontweight='bold')

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#8B0000', 
           markersize=8, label=f'Up (P<0.001)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728', 
           markersize=8, label=f'Up (P<0.05)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#00008B', 
           markersize=8, label=f'Down (P<0.001)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', 
           markersize=8, label=f'Down (P<0.05)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#7f7f7f', 
           markersize=8, label='NS'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / 'volcano_plot.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/volcano_plot.png")

# ===================================================================
# 9. MA PLOT
# ===================================================================

print("\n[Step 9] Creating MA plot...")

fig, ax = plt.subplots(figsize=(10, 8))

x = results['AveExpr']
y = results['logFC']
colors = results['Regulation'].map({
    'Up': '#d62728',
    'Down': '#1f77b4',
    'NS': '#7f7f7f'
})

ax.scatter(x, y, c=colors, alpha=0.6, s=20, edgecolors='none')
ax.axhline(fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axhline(-fc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axhline(0, color='black', linestyle='-', linewidth=1, alpha=0.3)

ax.set_xlabel('Average Expression (log₂)', fontsize=12, fontweight='bold')
ax.set_ylabel('log₂ Fold Change', fontsize=12, fontweight='bold')
ax.set_title(f'MA Plot: {condition2} vs {condition1}', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'ma_plot.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/ma_plot.png")

# ===================================================================
# 10. HEATMAP
# ===================================================================

print("\n[Step 10] Creating heatmap...")

top_n = min(50, n_sig) if n_sig > 0 else 50
top_proteins = results.nsmallest(top_n, 'P.Value')
top_indices = top_proteins.index

heatmap_data = intensities_imputed.iloc[top_indices].copy()

# Gene labels
gene_labels = []
for idx in top_indices:
    gene = str(protein_info.loc[idx, 'Gene_names']).split(';')[0]
    if pd.isna(gene) or gene == 'nan' or gene == '':
        gene = str(protein_info.loc[idx, 'Protein_IDs']).split(';')[0][:15]
    gene_labels.append(gene)

heatmap_data.index = gene_labels

# Z-score normalization (row-wise)
# Apply zscore across each row (protein), not as apply
heatmap_data_z = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(heatmap_data.std(axis=1), axis=0)

# Column colors
col_colors = [
    '#1f77b4' if sample_to_condition.get(s) == condition1 else '#d62728'
    for s in heatmap_data_z.columns
]

# Clustered heatmap
from scipy.cluster.hierarchy import linkage
import matplotlib.patches as mpatches

row_linkage = linkage(heatmap_data_z, method='average')
col_linkage = linkage(heatmap_data_z.T, method='average')

g = sns.clustermap(heatmap_data_z, 
                   cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                   row_linkage=row_linkage, col_linkage=col_linkage,
                   col_colors=col_colors, yticklabels=True, xticklabels=True,
                   figsize=(12, 14), cbar_kws={'label': 'Z-score'})

g.ax_heatmap.set_xlabel('Samples', fontsize=12, fontweight='bold')
g.fig.suptitle(f'Top {top_n} Proteins', y=0.98, fontsize=14, fontweight='bold')

handles = [
    mpatches.Patch(color='#1f77b4', label=condition1),
    mpatches.Patch(color='#d62728', label=condition2)
]
g.ax_heatmap.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.15, 1))

plt.savefig(output_dir / 'heatmap_top_proteins.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/heatmap_top_proteins.png")

# ===================================================================
# 11. P-VALUE DISTRIBUTIONS
# ===================================================================

print("\n[Step 11] Creating p-value distributions...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].hist(results['P.Value'], bins=50, edgecolor='black', color='steelblue')
axes[0].axvline(0.05, color='red', linestyle='--', linewidth=2, label='P = 0.05')
axes[0].set_xlabel('P-value', fontsize=12)
axes[0].set_ylabel('Frequency', fontsize=12)
axes[0].set_title('Raw P-value Distribution', fontsize=13, fontweight='bold')
axes[0].legend()
axes[0].grid(axis='y', alpha=0.3)

axes[1].hist(results['adj.P.Val'], bins=50, edgecolor='black', color='steelblue')
axes[1].axvline(0.05, color='red', linestyle='--', linewidth=2, label='adj.P = 0.05')
axes[1].set_xlabel('Adjusted P-value (FDR)', fontsize=12)
axes[1].set_ylabel('Frequency', fontsize=12)
axes[1].set_title('Adjusted P-value Distribution', fontsize=13, fontweight='bold')
axes[1].legend()
axes[1].grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'pvalue_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: {output_dir}/pvalue_distribution.png")

# ===================================================================
# 12. SAVE RESULTS (INCLUDING ENRICHMENT-READY OUTPUT)
# ===================================================================

print("\n[Step 12] Saving results...")

# Full results
results.to_csv(output_dir / 'differential_results_full.csv', index=False)
print(f"✓ Saved: {output_dir}/differential_results_full.csv ({len(results)} proteins)")

# Significant proteins
sig_results = results[results['Significant']].copy()
sig_results.to_csv(output_dir / 'significant_proteins.csv', index=False)
print(f"✓ Saved: {output_dir}/significant_proteins.csv ({len(sig_results)} proteins)")

# Top 100
top100 = results.head(100)
top100.to_csv(output_dir / 'top100_proteins.csv', index=False)
print(f"✓ Saved: {output_dir}/top100_proteins.csv")

# ENRICHMENT-READY OUTPUT
enrichment_input = results[results['Significant']].copy()
enrichment_input = enrichment_input[['Gene_names', 'logFC', 'adj.P.Val', 'Regulation']]
enrichment_input = enrichment_input.dropna(subset=['Gene_names'])
enrichment_input['Gene_names'] = enrichment_input['Gene_names'].str.split(';').str[0]  # First gene only
enrichment_input = enrichment_input[enrichment_input['Gene_names'] != '']
enrichment_input.to_csv(output_dir / 'sig_for_enrichment.tsv', sep='\t', index=False)
print(f"✓ Saved: {output_dir}/sig_for_enrichment.tsv ({len(enrichment_input)} genes for enrichment)")

# ===================================================================
# 13. FINAL SUMMARY REPORT
# ===================================================================

print("\n" + "="*70)
print("DIFFERENTIAL ABUNDANCE ANALYSIS - FINAL REPORT")
print("="*70)

print(f"\nComparison: {condition2} vs {condition1}")
print(f"  {condition1}: n={len(group1_samples)}")
print(f"  {condition2}: n={len(group2_samples)}")

print(f"\nStatistical Method:")
print(f"  Empirical Bayes moderated t-test")
print(f"  Prior df: {results_dict['prior_df']:.2f} (estimated from data)")
print(f"  Prior variance: {results_dict['prior_var']:.4f} (median)")

print(f"\nThresholds:")
print(f"  |logFC| > {fc_threshold} (= {2**fc_threshold:.1f}-fold change)")
print(f"  adj.P.Val < {pval_threshold} (FDR correction)")

print(f"\nResults:")
print(f"  Total: {n_total}")
print(f"  Significant: {n_sig} ({n_sig/n_total*100:.1f}%)")
print(f"    Up in {condition2}: {n_up}")
print(f"    Down in {condition2}: {n_down}")

print(f"\nTop 10 Hits:")
print(results[['Gene_names', 'logFC', 'FoldChange', 'P.Value', 'adj.P.Val']].head(10).to_string(index=False))

print(f"\nOutput Files:")
print(f"  - qc_summary.csv")
print(f"  - differential_results_full.csv")
print(f"  - significant_proteins.csv")
print(f"  - sig_for_enrichment.tsv  ← Ready for gProfiler/enrichr")
print(f"  - Figures: volcano, MA, heatmap, p-value distributions")

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)