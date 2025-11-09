#!/usr/bin/env python3
"""
Day 5: Target/Biomarker Prioritization Pipeline
Integrates all previous analyses to rank therapeutic targets and diagnostic biomarkers

Scoring Criteria:
1. Statistical significance (fold change + p-value)
2. Biological relevance (mitochondrial, pathways)
3. Network centrality (hub proteins)
4. Disease association (literature, databases)
5. Druggability (known drug targets)
6. Tissue specificity (brain/muscle expression)
"""
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

# Create output directory
output_dir = Path('biomarker_results')
output_dir.mkdir(exist_ok=True)

print("="*70)
print("DAY 5: TARGET/BIOMARKER PRIORITIZATION")
print("="*70)

# ==============================================================================
# 1. LOAD ALL PREVIOUS RESULTS
# ==============================================================================

print("\n[Step 1] Loading all previous analyses...\n")

# Differential results
diff_results = pd.read_csv('differential_results/differential_results_full.csv')
sig_proteins = pd.read_csv('differential_results/significant_proteins.csv')

print(f"✓ Differential: {len(diff_results)} proteins, {len(sig_proteins)} significant")

# Enrichment results
try:
    gprofiler_all = pd.read_csv('v3_enrichment_results/gprofiler_all_significant.csv')
    print(f"✓ Enrichment: {len(gprofiler_all)} pathways")
except:
    gprofiler_all = pd.DataFrame()
    print("⚠ No enrichment results found")

# Mitochondrial proteins
try:
    mito_proteins = pd.read_csv('v3_enrichment_results/mitochondrial_proteins.csv')
    print(f"✓ Mitochondrial: {len(mito_proteins)} proteins")
except:
    mito_proteins = pd.DataFrame()
    print("⚠ No mitochondrial results found")

# STRING network
try:
    string_stats = pd.read_csv('new_1_enrichment_results/string_network_statistics.csv')
    print(f"✓ STRING network: {len(string_stats)} proteins")
except:
    string_stats = pd.DataFrame()
    print("⚠ No STRING results found")

print()

# ==============================================================================
# 2. MULTI-CRITERIA SCORING SYSTEM
# ==============================================================================

print("[Step 2] Building multi-criteria scoring system...\n")

# Start with significant proteins
biomarker_candidates = sig_proteins.copy()

# Extract first gene name for scoring
biomarker_candidates['Primary_Gene'] = biomarker_candidates['Gene_names'].str.split(';').str[0]

print(f"Evaluating {len(biomarker_candidates)} candidate biomarkers\n")

# --- Score 1: Statistical Strength (0-20 points) ---
print("Calculating Score 1: Statistical Strength...")

# Fold change score (0-10 points)
fc_abs = np.abs(biomarker_candidates['logFC'])
biomarker_candidates['FC_Score'] = np.clip(
    (fc_abs - 1) / 2 * 10,  # Scale: 1-3 logFC → 0-10 points
    0, 10
)

# P-value score (0-10 points)
biomarker_candidates['Pval_Score'] = np.clip(
    -np.log10(biomarker_candidates['adj.P.Val']),
    0, 10
)

biomarker_candidates['Statistical_Score'] = (
    biomarker_candidates['FC_Score'] + biomarker_candidates['Pval_Score']
)

print(f"  Range: {biomarker_candidates['Statistical_Score'].min():.1f} - {biomarker_candidates['Statistical_Score'].max():.1f}")

# --- Score 2: Mitochondrial Relevance (0-15 points) ---
print("Calculating Score 2: Mitochondrial Relevance...")

# MitoCarta 3.0 gene list (download from broad -institute, expanded)
mitocarta_genes = [ 
    'CYC1', 'SDHB', 'COQ7', 'SDHA', 'UQCRC1', 'COQ5', 'PDHA1',
    'COQ9', 'MRPL12', 'ATP5F1D', 'COX5A', 'ISCA2', 'PMPCB',
    'UQCRFS1', 'ATP5F1A', 'OGDH', 'PDHB', 'UQCRC2', 'SDHD',
    'MRPS35', 'UQCRQ', 'MRPL53', 'DBT', 'PDK4', 'MDH2',
    'MRPS27', 'CS', 'GRPEL1', 'DLAT', 'LRPPRC', 'DLST', 'PDHX', 'GFM1',
    'MPC2', 'NDUFS1', 'MRPL46', 'ATP5F1E', 'SLC25A3', 'MRPS23',
    'FH', 'PMPCA', 'ATP5F1B', 'SDHAF4', 'UQCR10', 'ISCA1',
    'SUCLA2', 'COQ3', 'IARS2', 'MRPS15', 'IDH3A', 'COX11',
    'ETFDH', 'TIMM10', 'MRPL34', 'MRPL2', 'BCKDHA', 'UQCRH',
    'HIGD2A', 'ATP5PO', 'ECHS1', 'LETMD1', 'COX6A1', 'COX15',
    'AFG3L2', 'HADHA', 'ETFA', 'NDUFS7', 'CPT2', 'BCKDHB',
    'IDH3B', 'LARS2', 'ACADS', 'LETM1', 'ATP5ME', 'OPA1',
    'AUH', 'SUCLG1', 'NDUFV2', 'COQ6', 'MRPL43', 'ABHD11',
    'ATP5PF', 'NDUFB8', 'LONP1', 'DLD', 'AIFM1', 'ECHDC3',
    'APOOL', 'MRPL10', 'MRRF', 'NDUFS8', 'ACADM', 'IMMT',
    'TIMM9', 'SLC25A4', 'SAMM50', 'NDUFS2', 'NDUFV1', 'ACO2',
    'SUPV3L1', 'FECH', 'MTIF2', 'PHB', 'HIBCH', 'MRPS2', 'HSPA9', 'SURF1',
    'PRDX3', 'GHITM', 'GUF1', 'TIMM13', 'LYRM4', 'MRPL16', 'MRPL40', 'IDH3G',
    'SDHC', 'NDUFB5', 'SDHAF2', 'COQ10A', 'GATD3A',
    'TXN2', 'MRPS18A', 'COX6C', 'NDUFB9', 'MTCH2', 'NDUFA6', 'SLC25A20', 'MRPL1', 'TIMM44',
    'COX17', 'MICOS10', 'SOD2', 'COX6B1', 'VDAC1', 'CLPP', 'HADH', 'ACADL',
    'ACAT1', 'MALSU1', 'CLPX', 'NDUFS4', 'MRPL4', 'C1QBP', 'PITRM1',
    'UQCC1', 'MRPL55', 'MRM1', 'MECR', 'MRPL44', 'HSDL2', 'MRPS14',
    'ABCB8', 'ATPAF2', 'NDUFA9', 'COA6', 'POLRMT', 'ISCU', 'RTN4IP1',
    'OGDHL', 'ATP5F1C', 'BCKDK', 'GFM2', 'NDUFAF4', 'SLC25A35',
    'SLC25A11', 'TIMM8B', 'MCAT', 'IBA57', 'PHB2', 'DAP3',
    'CMC2', 'SDHAF1', 'MRPS17', 'ADHFE1', 'NFS1', 'NDUFAF5',
    'MRPL17', 'TACO1', 'MFN1', 'PPTC7', 'MRPL11', 'COX4I1',
    'NDUFS6', 'MRPL13', 'OXCT1', 'COX5B', 'PDK2', 'ALAS1',
    'MRPL33', 'MTO1', 'LIAS', 'NDUFA5', 'NDUFB6', 'MTX2',
    'SUCLG2', 'FDX1', 'SLC25A1', 'MRPS28', 'TIMM50', 'VDAC2',
    'SLC25A42', 'TOMM7', 'ECH1', 'GATB', 'BCS1L', 'ERAL1',
    'CMC1', 'MRPL28', 'TSFM', 'FXN', 'NFU1', 'YARS2',
    'MRPL21', 'ALDH4A1', 'TOMM40L', 'NDUFAF1', 'NDUFB10', 'ACADSB',
    'MRPS9', 'COQ8A', 'ACADVL', 'COQ4', 'MRPL15', 'NDUFA2',
    'MRPL24', 'HSPD1', 'NDUFS5', 'AK3', 'CYCS', 'MIPEP',
    'LYRM7', 'CRAT', 'PCCB', 'MRPS7', 'MRPL3', 'PRODH',
    'PCCA', 'MCCC1', 'MRPS12', 'CLPB', 'PDK1', 'MRPL49',
    'COX7A2', 'TMEM126A', 'ECHDC2', 'HCCS', 'HIBADH', 'MRPL19',
    'MRPL36', 'SLC25A30', 'BDH1', 'FARS2', 'ABCB7', 'MTX1',
    'NDUFA7', 'TIMM17A', 'ALDH9A1', 'MRPS18C', 'MARS2', 'ALDH6A1',
    'FDXR', 'GATC', 'TRAP1', 'ACAD8', 'ALDH2', 'PPIF',
    'TIMM22', 'IVD', 'L2HGDH', 'ETHE1', 'MRPL20', 'SLC25A5',
    'SLC25A12', 'MRPS21', 'TOMM22', 'ACAA2', 'MRPL30', 'DNAJA3',
    'NDUFB2', 'MRPS34', 'ETFB', 'AFG1L', 'ATP5IF1', 'ATP5PB',
    'COX7C', 'CHCHD3', 'COQ10B', 'ACSF3', 'POLDIP2', 'SLC25A10',
    'SLC25A13', 'PDK3', 'ME3', 'MRPL22', 'IDH2', 'GCDH',
    'MRPL47', 'PPA2', 'MRPL9', 'CHCHD10', 'WARS2', 'SLC25A19',
    'CBR4', 'SMDT1', 'HAGH', 'COX7A1', 'MTG1', 'COX14',
    'NDUFA12', 'MRPS16', 'MTERF2', 'ETFRF1', 'GPD2', 'NDUFA10',
    'NDUFC2', 'SLIRP', 'MT-ATP6', 'MT-CO2', 'TIMM21', 'NDUFA8',
    'GLS', 'NDUFB11', 'DHTKD1', 'FAHD2A', 'HSD17B8', 'HINT2',
    'MRPS5', 'C6orf136', 'SPRYD4', 'LIPT2', 'DECR1', 'MRPS26',
    'SLC25A15', 'NDUFV3', 'BPHL', 'STOML2', 'MSRB2', 'LACTB',
    'TOMM40', 'SLC25A25', 'ADCK5', 'ACAD10', 'VWA8', 'CCDC90B',
    'DNAJC11', 'ALKBH7', 'GLDC', 'MPC1L', 'THEM4', 'SSBP1',
    'MRPL27', 'HSCB', 'MRPS10', 'AK4', 'ATAD3A', 'MRPL37',
    'PTGES2', 'TXNRD2', 'ACSM5', 'CISD1', 'MRPL41', 'MTARC2',
    'MRPS6', 'MRPL23', 'SQOR', 'SCO1', 'GCAT', 'MTHFD1L',
    'ECI2', 'UQCRB', 'ATP5PD', 'ATP5MG', 'CCDC58', 'MCCC2',
    'MCEE', 'MOCS1', 'MMUT', 'SLC25A51', 'SPG7', 'ATP5MC3',
    'ATAD1', 'OXNAD1', 'MTG2', 'DARS2', 'MRPL51', 'PC',
    'CLYBL', 'PNPLA8', 'MLYCD', 'PYURF', 'TIMM17B', 'MRPL32',
    'SLC25A26', 'PRODH2', 'NIT2', 'FDX2', 'ME2', 'TUFM',
    'COX6A2', 'SLC30A9', 'TAMM41', 'TIMM23', 'RMND1', 'MFN2',
    'DNLZ', 'CHCHD4', 'NDUFAB1', 'ACSM1', 'COA5', 'COQ2',
    'PDSS2', 'EHHADH', 'NDUFS3', 'ABCB10', 'CKMT2', 'MRPS25',
    'ACAD9', 'COX16', 'FAM210A', 'ACOT13', 'CPT1A', 'DHRS4',
    'PRELID2', 'CARS2', 'NIPSNAP2', 'NDUFB7', 'MICOS13', 'GLRX2',
    'MRPS30', 'SCO2', 'HSD17B10', 'NARS2', 'RPUSD3', 'SUOX',
    'SARDH', 'RIDA', 'COX8A', 'SPR', 'MPV17', 'ZADH2',
    'MRPL18', 'ALDH1B1', 'CHCHD1', 'TIMM8A', 'VDAC3', 'SLC25A24',
    'MPC1', 'SLC25A22', 'SLC25A31', 'SLC25A18', 'TRNT1', 'IDE',
    'GSTK1', 'VARS2', 'HSPE1', 'MRPL54', 'MRPL38', 'MRPL58',
    'MTCH1', 'TRIAP1', 'PTCD3', 'COX10', 'AKAP1', 'MRM3',
    'AK2', 'TMEM70', 'NOA1', 'MCUR1', 'TMEM14C', 'NDUFC1',
    'MRPL35', 'NDUFA3', 'PDSS1', 'SLC25A16', 'FMC1', 'ACSM3',
    'NDUFAF6', 'UCP1', 'SFXN1', 'GLRX5', 'PTRH1', 'MRPS11',
    'XPNPEP3', 'PDP1', 'ECI1', 'SLC25A21', 'SLC25A29', 'SLC25A45',
    'SDHAF3', 'PRELID3B', 'SLC25A39', 'NDUFB3', 'REXO2', 'GATM',
    'TIMM10B', 'SARS2', 'CHDH', 'MRPL57', 'C8orf82', 'EARS2',
    'GRHPR', 'MTFMT', 'PTCD2', 'ECHDC1', 'HARS2', 'NDUFA4',
    'NDUFA13', 'NAXE', 'TMEM11', 'SELENOO', 'NAGS', 'HSD17B4',
    'NDUFAF7', 'OAT', 'MRPL50', 'COX7B', 'RARS2', 'OXLD1',
    'LACTB2', 'NDUFA11', 'AURKAIP1', 'CHCHD7', 'SHMT2', 'MTHFD2',
    'RMDN1', 'CPS1', 'MRPL14', 'MTRES1', 'RDH13', 'TRMT2B',
    'MRPS24', 'COA3', 'GSR', 'ENDOG', 'BCAT2', 'TEFM',
    'ACSF2', 'TRMU', 'MRPS22', 'LAP3', 'GLUD1', 'HMGCL',
    'FHIT', 'MRPS33', 'NUBPL', 'LYRM2', 'MPV17L', 'GADD45GIP1',
    'TFAM', 'ALDH18A1', 'OXSM', 'KYAT3', 'TARS2', 'TBRG4',
    'TMEM65', 'ROMO1', 'PRDX5', 'DNAJC15', 'MAIP1', 'HIGD1A',
    'C16orf91', 'RHOT2', 'SLC25A44', 'HAO2', 'OXA1L', 'SLC25A46',
    'RDH14', 'ATP5MPL', 'CA5B', 'COASY', 'CISD3', 'GSTZ1',
    'COX6B2', 'GLOD4', 'ACAA1', 'CROT', 'PDF', 'COX19',
    'AMT', 'PDPR', 'NLN', 'BOLA1', 'ALAS2', 'GRPEL2',
    'ALDH7A1', 'HDHD5', 'GFER', 'IMMP2L', 'COX4I2', 'MRPL39',
    'AADAT', 'FUNDC2', 'LYPLA1', 'PDHA2', 'MTRF1L', 'GTPBP3',
    'MTHFD2L', 'YME1L1', 'ACOT2', 'COA8', 'MMAB', 'PPOX',
    'ABCD3', 'SIRT3', 'LDHD', 'MRPL45', 'C12orf65', 'IMMP1L',
    'POLG', 'ATPAF1', 'BLOC1S1', 'ALDH1L2', 'COX18', 'CPOX',
    'NNT', 'NDUFAF8', 'ATP23', 'DCAKD', 'ATP5MF', 'SFXN3',
    'AASS', 'DIABLO', 'COX8C', 'PXMP2', 'AGXT2', 'GLS2',    
    'APOO', 'ACAD11', 'C5orf63', 'HOGA1', 'TTC19', 'NGRN', 'TCAIM',
    'GCSH', 'SIRT5', 'GPT2', 'TOMM6', 'FAM162A', 'ABCB6',
    'MACROD1', 'TST', 'MT-CO1', 'MTFR1L', 'SDR39U1', 'DMGDH',
    'SLC25A14', 'PRELID3A', 'CHCHD5', 'TOMM70', 'NIT1', 'FAHD1',
    'MCU', 'CYP27A1', 'MT-ND2', 'MT-ND4', 'MT-ND5', 'DELE1',
    'FASTKD1', 'MRPS31', 'ALDH1L1', 'PTCD1', 'HADHB', 'TRIT1',
    'CRLS1', 'GLYCTK', 'SLC25A32', 'AARS2', 'EXOG', 'CHCHD2',
    'KMO', 'SFXN5', 'ELAC2', 'GPX1', 'CPT1B', 'FIS1',
    'PRELID1', 'UCP3', 'SLC25A27', 'SLC25A40', 'SLC25A23', 'HEMK1',
    'SLC25A36', 'OMA1', 'LETM2', 'DHRS1', 'STOM', 'METTL5',
    'MSRA', 'NUDT8', 'NMNAT3', 'SLC25A38', 'MPV17L2', 'ACLY',
    'ABCD2', 'FLAD1', 'LIPT1', 'DGLUCY', 'PCK2', 'NADK2',
    'MPST', 'NUDT13', 'ACACA', 'ABAT', 'HMGCS2', 'ALDH5A1',
    'CKMT1A', 'PISD', 'NDUFB4', 'UQCR11', 'FAM136A', 'TOMM5',
    'HTRA2', 'NDUFA1', 'MT-CYB', 'GUK1', 'D2HGDH', 'OPA3',
    'AKR1B10', 'PHYH', 'NDUFAF2', 'ECSIT', 'PARL', 'MRPL42',
    'DNAJC4', 'PAM16', 'ALDH3A2', 'ACSL6', 'MRS2', 'ABHD10',
    'HINT3', 'PXMP4', 'OSGEPL1', 'NUDT2', 'CYB5B', 'GRSF1',
    'MRPS18B', 'DNAJC30', 'CAT', 'DGUOK', 'ACOT9', 'ACSS1',
    'RFK', 'STARD7', 'CHPT1', 'MT-ATP8', 'TIMMDC1', 'SPATA19',
    'HTATIP2', 'AKR7A2', 'DUT', 'MT-CO3', 'DHRS7B', 'MTRF1',
    'TMEM143', 'NT5M', 'RHOT1', 'NME4', 'MTHFS', 'TIMM29',
    'PNPO', 'MICU2', 'PCBD2', 'QRSL1', 'ADCK1', 'PRDX6',
    'NIPSNAP1', 'GARS1', 'ACSM2A', 'SLC25A33', 'LYRM9', 'RAB24',
    'PARS2', 'MRPL52', 'PUS1', 'GOT2', 'NAXD', 'MT-ND3',
    'ACCS', 'PTRH2', 'PGS1', 'OXR1', 'LYPLAL1', 'DNAJC28',
    'FPGS', 'TMEM205', 'SLC25A53', 'AIFM3', 'MGST1', 'MTFP1',
    'COA4', 'NIF3L1', 'PUSL1', 'SDSL', 'BAX', 'MAVS',
    'METTL17', 'ACACB', 'ARG2', 'MTPAP', 'CCDC51', 'SUGCT',
    'DNM1L', 'NRDC', 'TMEM126B', 'NSUN4', 'ACSS3', 'FAM210B',
    'RBFA', 'METAP1D', 'AHCYL1', 'NIPSNAP3A', 'MTARC1', 'FUNDC1',
    'ACSL1', 'COX20', 'HSDL1', 'CBR3', 'MTERF4', 'NUDT19',
    'HEBP1', 'PTPMT1', 'PGAM5', 'DTYMK', 'KARS1', 'STAR',
    'DNAJC19', 'OTC', 'COX7A2L', 'MTERF3', 'SFXN2', 'FASN',
    'PINK1', 'TMEM177', 'TXNRD1', 'THNSL1', 'FDPS', 'IDI1',
    'TOP3A', 'FTMT', 'LDHAL6B', 'MTFR1', 'BCL2L13', 'TRMT10C',
    'ACOT7', 'TFB1M', 'ARF5', 'MCUB', 'TWNK', 'MMADHC',
    'PYCR2', 'TSTD1', 'LYRM1', 'OSBPL1A', 'RECQL4', 'PET100',
    'GDAP1', 'MTIF3', 'UNG', 'GLYAT', 'ISOC2', 'ATP5MC1',
    'QDPR', 'COQ8B', 'GTPBP10', 'TSPO', 'FTH1', 'TRMT1',
    'ADCK2', 'TSTD3', 'MRPS36', 'PDP2', 'AGPAT5', 'MCRIP2',
    'UCP2', 'SLC25A37', 'SLC25A28', 'SLC25A48', 'SLC25A43', 'SLC25A47',
    'SLC25A41', 'SLC25A34', 'EPHX2', 'NUDT5', 'MICU1', 'SOD1',
    'SERHL2', 'NDUFAF3', 'BOLA3', 'OCIAD2', 'RSAD1', 'CHCHD6',
    'FASTKD2', 'PEX11B', 'TMEM186', 'PIF1', 'FABP1', 'ABCD1',
    'PLPBP', 'DBI', 'PNPT1', 'TOP1MT', 'RNASEH1', 'ABCA9',
    'SCP2', 'NT5DC3', 'CMC4', 'LDHB', 'YBEY', 'HINT1',
    'MUL1', 'TK2', 'AGMAT', 'CMPK2', 'BID', 'ARL2',
    'SYNJ2BP', 'MFF', 'DCXR', 'NUDT9', 'AMACR', 'VPS13D',
    'NUDT6', 'DMAC2', 'MRM2', 'BAK1', 'DHODH', 'MAOB',
    'MGST3', 'ACP6', 'MRPL48', 'AGXT', 'CYB5R3', 'DDX28',
    'PNKD', 'MICU3', 'BNIP3', 'MTX3', 'MUTYH', 'GTPBP6',
    'CYP11A1', 'OCIAD1', 'METTL8', 'SMIM20', 'RPUSD4', 'SMIM8',
    'NME6', 'FKBP8', 'ACSM4', 'NLRX1', 'FOXRED1', 'TFB2M',
    'PPM1K', 'ATP5MC2', 'MT-ND1', 'TAZ', 'MSRB3', 'ACOT11',
    'NIPSNAP3B', 'NTHL1', 'MARCHF5', 'SFXN4', 'TRMT5', 'SLC25A6',
    'SLC25A52', 'OGG1', 'ACOD1', 'QTRT1', 'AIFM2', 'PRDX2',
    'FBXL4', 'GPAM', 'DMAC2L', 'RMDN3', 'PRXL2A', 'ATP5MD',
    'PET117', 'MAOA', 'C15orf61', 'POLQ', 'BCO2', 'IFI27',
    'CRYZ', 'PYCR1', 'TMLHE', 'TOMM20', 'YRDC', 'COA7',
    'DHRS2', 'MMAA', 'MIEF1', 'C15orf48', 'LIG3', 'BNIP3L',
    'FASTKD5', 'PARK7', 'COA1', 'COX7B2', 'MT-ND4L', 'SIRT4',
    'UQCC2', 'USP30', 'PDE12', 'NME3', 'DHX30', 'FASTK',
    'METTL15', 'CYP24A1', 'NOCT', 'PANK2', 'CA5A', 'NT5DC2',
    'ATAD3B', 'TRUB2', 'PRKACA', 'COMTD1', 'CPT1C', 'C3orf33',
    'MIGA2', 'SND1', 'GLUD2', 'HPDL', 'ANGEL2', 'ANTKMT',
    'UQCC3', 'PRDX4', 'CDK5RAP1', 'PAICS', 'MGARP', 'HDHD3',
    'CYP11B2', 'AGPAT4', 'RCC1L', 'FAM185A', 'PMAIP1', 'TDRKH',
    'SLC8B1', 'CYP27B1', 'POLG2', 'THEM5', 'BIK', 'MTERF1',
    'MTFR2', 'GPX4', 'TOMM20L', 'COMT', 'MYG1', 'SPATA20',
    'OXCT2', 'CCDC127', 'CRY1', 'ACSM2B', 'MYO19', 'MGME1',
    'SEPTIN4', 'EFHD1', 'RAB5IF', 'APEX1', 'SPTLC2', 'DMAC1',
    'ATPSCKMT', 'NEU4', 'CASP9', 'AGK', 'MCL1', 'DUS2',
    'PLD6', 'TOMM34', 'PRORP', 'PLSCR3', 'CKMT1B', 'NDUFB1',
    'ATP5MF-PTCD1', 'CASP8', 'BCL2L1', 'EXD2', 'ARMC10', 'MIEF2',
    'THG1L', 'GOLPH3', 'PREPL', 'FASTKD3', 'CEP89', 'GPAT2',
    'NBR1', 'BOK', 'FKBP10', 'SETD9', 'CYP11B1', 'CSKMT',
    'MIGA1', 'RPIA', 'BBC3', 'BAD', 'BCL2', 'BCL2A1',
    'BCL2L2', 'CASP3', 'DMPK', 'DNA2', 'MT-ND6', 'PRKN',
    'ALKBH1', 'PICK1', 'BCL2L10', 'AKAP10', 'STYXL1', 'NSUN2',
    'ADCY10', 'SPHK2', 'NSUN3', 'PABPC5', 'PRSS35', 'PRIMPOL',
    'SERAC1', 'STX17', 'TRMT61B', 'PDE2A', 'PLGRKT', 'ETFBKMT',
    'RTL10', 'BCL2L11', 'SPIRE1', 'C2orf69', 'ARMCX2', 'ARMCX1',
    'POLB', 'SNAP29', 'ARMCX3', 'ARMCX6', 'METTL4', 'SPHKAP',
    'NAT8L', 'MCCD1', 'PIGBOS1', 'HTD2', 'RP11_469A15.2', 'SURF1'
]

# Check if gene is mitochondrial
biomarker_candidates['Is_Mitochondrial'] = biomarker_candidates['Primary_Gene'].isin(mitocarta_genes)

# Respiratory chain complexes (higher priority)
complex_I = [g for g in mitocarta_genes if g.startswith(('NDUF', 'MT-ND'))]
complex_II = [g for g in mitocarta_genes if g.startswith('SDH')]
complex_III = [g for g in mitocarta_genes if g.startswith(('UQCR', 'CYC'))]
complex_IV = [g for g in mitocarta_genes if g.startswith(('COX', 'MT-CO'))]
complex_V = [g for g in mitocarta_genes if g.startswith(('ATP5', 'MT-ATP'))]

def get_mito_score(gene):
    if gene in complex_I:
        return 15  # Highest priority (Complex I most affected in Leigh)
    elif gene in complex_IV:
        return 12  # High priority
    elif gene in complex_V:
        return 10
    elif gene in complex_II or gene in complex_III:
        return 8
    elif gene in mitocarta_genes:
        return 5  # General mitochondrial
    else:
        return 0

biomarker_candidates['Mito_Score'] = biomarker_candidates['Primary_Gene'].apply(get_mito_score)

print(f"  Mitochondrial proteins: {biomarker_candidates['Is_Mitochondrial'].sum()}")
print(f"  Complex I genes: {sum(biomarker_candidates['Primary_Gene'].isin(complex_I))}")

# --- Score 3: Network Centrality (0-15 points) ---
print("Calculating Score 3: Network Centrality...")

if not string_stats.empty:
    # Normalize degree to 0-15 scale
    max_degree = string_stats['Degree'].max()
    if max_degree > 0:
        degree_dict = dict(zip(string_stats['Gene'], string_stats['Degree']))
        biomarker_candidates['Network_Score'] = biomarker_candidates['Primary_Gene'].map(degree_dict).fillna(0)
        biomarker_candidates['Network_Score'] = (
            biomarker_candidates['Network_Score'] / max_degree * 15
        )
    else:
        biomarker_candidates['Network_Score'] = 0
    
    n_networked = (biomarker_candidates['Network_Score'] > 0).sum()
    print(f"  Proteins in network: {n_networked}")
else:
    biomarker_candidates['Network_Score'] = 0
    print("  No network data available")
# --- Score 4: Pathway Enrichment (WORKING VERSION) ---
print("Calculating Score 4: Pathway Enrichment...")

if not gprofiler_all.empty:
    # Score based on mitochondrial relevance of pathway terms
    mitochondrial_terms = {
        'high': ['mitochondri', 'respiratory chain', 'electron transport', 'oxidative phosphorylation',
                'complex i', 'complex ii', 'complex iii', 'complex iv', 'complex v',
                'atp synthesis', 'atp metabolic process', 'energy derivation'],
        'medium': ['tca cycle', 'krebs cycle', 'citrate cycle', 'fatty acid oxidation',
                  'beta oxidation', 'reactive oxygen species', 'redox process',
                  'cellular respiration', 'aerobic respiration'],
        'low': ['metabolic process', 'catabolic process', 'biosynthetic process',
               'cellular homeostasis', 'ion transport', 'transmembrane transport']
    }
    
    # Calculate pathway relevance score for each term
    pathway_scores = {}
    
    for _, row in gprofiler_all.iterrows():
        term_name = str(row['Term_Name']).lower()
        pval = row['P_value']
        source = row['Source']
        
        # Base score from term relevance
        base_score = 0
        if any(term in term_name for term in mitochondrial_terms['high']):
            base_score = 12
        elif any(term in term_name for term in mitochondrial_terms['medium']):
            base_score = 8
        elif any(term in term_name for term in mitochondrial_terms['low']):
            base_score = 4
        
        # Adjust for significance
        if base_score > 0:
            if pval < 0.001:
                base_score += 3
            elif pval < 0.01:
                base_score += 2
            elif pval < 0.05:
                base_score += 1
            
            # Database preference
            if source in ['KEGG', 'REAC']:
                base_score += 2
            elif source.startswith('GO:'):
                base_score += 1
            
            base_score = min(base_score, 15)
            
            # Store the maximum score for this pathway
            pathway_scores[term_name] = max(pathway_scores.get(term_name, 0), base_score)
    
    # Assign pathway scores to proteins based on mitochondrial status
    # Mitochondrial proteins get high pathway scores, others get medium if significant
    def assign_pathway_score(row):
        if row['Is_Mitochondrial']:
            # Mitochondrial proteins get high pathway scores
            if row['Mito_Score'] >= 10:  # Complex I/IV proteins
                return min(12 + row['Statistical_Score'] / 20 * 3, 15)
            else:  # Other mitochondrial proteins
                return min(8 + row['Statistical_Score'] / 20 * 4, 12)
        else:
            # Non-mitochondrial proteins get lower scores based on statistical significance
            return min(row['Statistical_Score'] / 20 * 8, 8)
    
    biomarker_candidates['Pathway_Score'] = biomarker_candidates.apply(assign_pathway_score, axis=1)
    
    # Show pathway term analysis
    if pathway_scores:
        print(f"  Found {len(pathway_scores)} relevant pathway terms")
        high_score_terms = {k: v for k, v in pathway_scores.items() if v >= 10}
        if high_score_terms:
            print("  High-scoring mitochondrial pathways:")
            for term, score in list(high_score_terms.items())[:3]:
                print(f"    • {term[:50]}... ({score} pts)")
    
    n_high_pathway = (biomarker_candidates['Pathway_Score'] >= 10).sum()
    n_medium_pathway = ((biomarker_candidates['Pathway_Score'] >= 5) & 
                        (biomarker_candidates['Pathway_Score'] < 10)).sum()
    
    print(f"  Proteins with high pathway scores (≥10): {n_high_pathway}")
    print(f"  Proteins with medium pathway scores (5-9): {n_medium_pathway}")
    print(f"  Max pathway score: {biomarker_candidates['Pathway_Score'].max():.1f}")
    
else:
    biomarker_candidates['Pathway_Score'] = 0
    print("  No pathway data available")

# --- Score 5: Known Disease Genes (0-20 points) ---
print("Calculating Score 5: Known Disease Association...")

# Known Leigh syndrome genes (from literature)
leigh_genes = [
    'SURF1', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS7', 'NDUFS8',
    'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA2', 'NDUFA9', 'NDUFA10', 'NDUFA11',
    'NDUFAF1', 'NDUFAF2', 'NDUFAF3', 'NDUFAF4', 'NDUFAF5', 'NDUFAF6',
    'SDHA', 'SDHAF1', 'COX10', 'COX15', 'SCO1', 'SCO2', 'LRPPRC',
    'PDHA1', 'PDHX', 'DLAT', 'DLD', 'MT-ATP6', 'MT-ND1', 'MT-ND2',
    'MT-ND3', 'MT-ND4', 'MT-ND5', 'MT-ND6'
]

# Mitochondrial disease genes (broader)
mito_disease_genes = leigh_genes + [
    'POLG', 'TWNK', 'TYMP', 'TK2', 'DGUOK', 'MPV17', 'MGME1',
    'OPA1', 'MFN2', 'AFG3L2', 'SPG7', 'ATAD3A'
]

def disease_score(gene):
    if gene in leigh_genes:
        return 20  # Known Leigh gene
    elif gene in mito_disease_genes:
        return 10  # Other mitochondrial disease
    else:
        return 0

biomarker_candidates['Disease_Score'] = biomarker_candidates['Primary_Gene'].apply(disease_score)

n_known = (biomarker_candidates['Disease_Score'] > 0).sum()
print(f"  Known disease genes: {n_known}")

# --- Score 6: Druggability (0-15 points) ---
print("Calculating Score 6: Druggability...")

# Druggable protein classes
druggable_families = {
    'kinases': ['kinase', 'MAPK', 'CDK', 'PKC', 'AKT', 'ERK'],
    'receptors': ['receptor', 'GPCR', 'EGFR', 'VEGFR'],
    'enzymes': ['synthase', 'ase', 'oxidase', 'reductase', 'dehydrogenase'],
    'transporters': ['SLC', 'ABC', 'transporter', 'channel'],
    'transcription_factors': ['transcription factor', 'TF', 'nuclear receptor']
}

def druggability_score(protein_name, gene_name):
    score = 0
    name_lower = str(protein_name).lower() + ' ' + str(gene_name).lower()
    
    # Check protein families
    for family, keywords in druggable_families.items():
        if any(kw.lower() in name_lower for kw in keywords):
            score += 5
            break
    
    # Membrane/secreted proteins (easier to target)
    if any(word in name_lower for word in ['membrane', 'receptor', 'channel', 'transporter']):
        score += 3
    
    # Enzymatic activity (druggable)
    if any(word in name_lower for word in ['enzyme', 'catalytic', 'ase activity']):
        score += 2
    
    return min(score, 15)

biomarker_candidates['Druggability_Score'] = biomarker_candidates.apply(
    lambda x: druggability_score(x['Protein_names'], x['Primary_Gene']),
    axis=1
)

n_druggable = (biomarker_candidates['Druggability_Score'] > 0).sum()
print(f"  Potentially druggable: {n_druggable}")

print()

# ==============================================================================
# 3. CALCULATE TOTAL SCORE & RANK
# ==============================================================================

print("[Step 3] Calculating total scores and ranking...\n")

# Total score (max 100 points)
biomarker_candidates['Total_Score'] = (
    biomarker_candidates['Statistical_Score'] +
    biomarker_candidates['Mito_Score'] +
    biomarker_candidates['Network_Score'] +
    biomarker_candidates['Pathway_Score'] +
    biomarker_candidates['Disease_Score'] +
    biomarker_candidates['Druggability_Score']
)

# Rank biomarkers
biomarker_candidates = biomarker_candidates.sort_values('Total_Score', ascending=False)
biomarker_candidates['Rank'] = range(1, len(biomarker_candidates) + 1)

# Classification
def classify_biomarker(row):
    if row['Total_Score'] >= 70:
        return 'Tier 1: High Priority'
    elif row['Total_Score'] >= 50:
        return 'Tier 2: Medium Priority'
    elif row['Total_Score'] >= 30:
        return 'Tier 3: Low Priority'
    else:
        return 'Tier 4: Exploratory'

biomarker_candidates['Priority_Tier'] = biomarker_candidates.apply(classify_biomarker, axis=1)

# Summary
tier_counts = biomarker_candidates['Priority_Tier'].value_counts()
print("Priority Distribution:")
for tier, count in tier_counts.items():
    print(f"  {tier}: {count}")

print(f"\nTop 10 Biomarker Candidates:")
top10 = biomarker_candidates.head(10)
for idx, row in top10.iterrows():
    print(f"  {row['Rank']:2d}. {row['Primary_Gene']:15s} (Score: {row['Total_Score']:.1f}, {row['Priority_Tier']})")

print()

# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

print("[Step 4] Saving prioritization results...\n")

# Full ranking
output_cols = [
    'Rank', 'Primary_Gene', 'Gene_names', 'Protein_names',
    'logFC', 'FoldChange', 'adj.P.Val', 'Regulation',
    'Total_Score', 'Priority_Tier',
    'Statistical_Score', 'Mito_Score', 'Network_Score',
    'Pathway_Score', 'Disease_Score', 'Druggability_Score',
    'Is_Mitochondrial'
]

biomarker_ranking = biomarker_candidates[output_cols].copy()
biomarker_ranking.to_csv(output_dir / 'biomarker_ranking_full.csv', index=False)
print(f"✓ Saved: {output_dir}/biomarker_ranking_full.csv ({len(biomarker_ranking)} candidates)")

# Top 50 candidates
top50 = biomarker_ranking.head(50)
top50.to_csv(output_dir / 'biomarker_top50.csv', index=False)
print(f"✓ Saved: {output_dir}/biomarker_top50.csv")

# Tier 1 only (high priority)
tier1 = biomarker_ranking[biomarker_ranking['Priority_Tier'] == 'Tier 1: High Priority']
if len(tier1) > 0:
    tier1.to_csv(output_dir / 'biomarker_tier1_high_priority.csv', index=False)
    print(f"✓ Saved: {output_dir}/biomarker_tier1_high_priority.csv ({len(tier1)} targets)")

# Known Leigh genes found
known_leigh_found = biomarker_ranking[biomarker_ranking['Disease_Score'] == 20]
if len(known_leigh_found) > 0:
    known_leigh_found.to_csv(output_dir / 'known_leigh_genes_validated.csv', index=False)
    print(f"✓ Saved: {output_dir}/known_leigh_genes_validated.csv ({len(known_leigh_found)} genes)")

print()

# ==============================================================================
# 5. VISUALIZATION
# ==============================================================================

print("[Step 5] Creating prioritization visualizations...\n")

# --- Plot 1: Score Distribution Radar ---
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Top 10 radar chart
from math import pi

categories = ['Statistical', 'Mitochondrial', 'Network', 'Pathway', 'Disease', 'Druggability']
N = len(categories)

angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

ax = plt.subplot(221, polar=True)

colors = plt.cm.Set3(range(10))
for i, (idx, row) in enumerate(top10.head(5).iterrows()):
    values = [
        row['Statistical_Score'], row['Mito_Score'], row['Network_Score'],
        row['Pathway_Score'], row['Disease_Score'], row['Druggability_Score']
    ]
    values += values[:1]
    
    ax.plot(angles, values, 'o-', linewidth=2, label=row['Primary_Gene'], color=colors[i])
    ax.fill(angles, values, alpha=0.15, color=colors[i])

ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, size=9)
ax.set_ylim(0, 20)
ax.set_title('Top 5 Candidates: Multi-Criteria Scores', size=12, fontweight='bold', pad=20)
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=9)
ax.grid(True)

# Plot 2: Score component breakdown
ax2 = plt.subplot(222)

score_components = biomarker_ranking.head(20)[[
    'Primary_Gene', 'Statistical_Score', 'Mito_Score', 'Network_Score',
    'Pathway_Score', 'Disease_Score', 'Druggability_Score'
]]

score_components.set_index('Primary_Gene').plot(
    kind='barh', stacked=True, ax=ax2, 
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
)

ax2.set_xlabel('Total Score', fontweight='bold')
ax2.set_title('Top 20: Score Component Breakdown', fontweight='bold')
ax2.legend(title='Score Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
ax2.invert_yaxis()

# Plot 3: Total score distribution
ax3 = plt.subplot(223)

ax3.hist(biomarker_ranking['Total_Score'], bins=30, edgecolor='black', color='steelblue', alpha=0.7)
ax3.axvline(70, color='red', linestyle='--', linewidth=2, label='Tier 1 threshold')
ax3.axvline(50, color='orange', linestyle='--', linewidth=2, label='Tier 2 threshold')
ax3.axvline(30, color='yellow', linestyle='--', linewidth=2, label='Tier 3 threshold')

ax3.set_xlabel('Total Score', fontweight='bold')
ax3.set_ylabel('Number of Candidates', fontweight='bold')
ax3.set_title('Score Distribution Across All Candidates', fontweight='bold')
ax3.legend()
ax3.grid(axis='y', alpha=0.3)

# Plot 4: Priority tiers
ax4 = plt.subplot(224)

tier_colors = {
    'Tier 1: High Priority': '#2ecc71',
    'Tier 2: Medium Priority': '#f39c12',
    'Tier 3: Low Priority': '#3498db',
    'Tier 4: Exploratory': '#95a5a6'
}

tier_data = biomarker_ranking['Priority_Tier'].value_counts()
colors_ordered = [tier_colors[tier] for tier in tier_data.index]

wedges, texts, autotexts = ax4.pie(
    tier_data.values,
    labels=tier_data.index,
    autopct='%1.1f%%',
    colors=colors_ordered,
    startangle=90
)

for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontweight('bold')

ax4.set_title('Biomarker Priority Distribution', fontweight='bold')

plt.tight_layout()
plt.savefig(output_dir / 'biomarker_prioritization_summary.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Saved: {output_dir}/biomarker_prioritization_summary.png")

# --- Plot 2: Heatmap of Top 30 ---
fig, ax = plt.subplots(figsize=(14, 10))

top30_scores = biomarker_ranking.head(30)[[
    'Primary_Gene', 'Statistical_Score', 'Mito_Score', 'Network_Score',
    'Pathway_Score', 'Disease_Score', 'Druggability_Score'
]]

top30_scores_matrix = top30_scores.set_index('Primary_Gene')

# Normalize to 0-1 for heatmap
top30_normalized = top30_scores_matrix.div(top30_scores_matrix.max(axis=0), axis=1)

sns.heatmap(
    top30_normalized.T,
    cmap='RdYlGn',
    linewidths=0.5,
    cbar_kws={'label': 'Normalized Score'},
    yticklabels=['Statistical', 'Mitochondrial', 'Network', 'Pathway', 'Disease', 'Druggability'],
    ax=ax
)

ax.set_title('Top 30 Biomarker Candidates: Score Heatmap', fontsize=14, fontweight='bold', pad=15)
ax.set_xlabel('Candidate Gene', fontsize=12, fontweight='bold')
ax.set_ylabel('Score Category', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(output_dir / 'biomarker_heatmap_top30.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Saved: {output_dir}/biomarker_heatmap_top30.png")

# --- Plot 3: Volcano plot with priority tiers ---
fig, ax = plt.subplots(figsize=(12, 9))

tier_colors_scatter = {
    'Tier 1: High Priority': '#2ecc71',
    'Tier 2: Medium Priority': '#f39c12',
    'Tier 3: Low Priority': '#3498db',
    'Tier 4: Exploratory': '#95a5a6'
}

for tier in biomarker_ranking['Priority_Tier'].unique():
    tier_data = biomarker_ranking[biomarker_ranking['Priority_Tier'] == tier]
    ax.scatter(
        tier_data['logFC'],
        -np.log10(tier_data['adj.P.Val']),
        c=tier_colors_scatter[tier],
        label=tier,
        s=100,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5
    )

# Label top 10
for _, row in biomarker_ranking.head(10).iterrows():
    ax.annotate(
        row['Primary_Gene'],
        (row['logFC'], -np.log10(row['adj.P.Val'])),
        fontsize=9,
        fontweight='bold',
        xytext=(5, 5),
        textcoords='offset points',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5)
    )

ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axvline(1, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.axvline(-1, color='black', linestyle='--', linewidth=1, alpha=0.5)

ax.set_xlabel('log₂ Fold Change', fontsize=13, fontweight='bold')
ax.set_ylabel('-log₁₀(Adjusted P-value)', fontsize=13, fontweight='bold')
ax.set_title('Biomarker Candidates: Prioritization Tiers', fontsize=14, fontweight='bold')
ax.legend(title='Priority Tier', loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'biomarker_volcano_prioritized.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Saved: {output_dir}/biomarker_volcano_prioritized.png")

print()

# ==============================================================================
# 6. GENERATE COMPREHENSIVE REPORT
# ==============================================================================

print("[Step 6] Generating comprehensive biomarker report...\n")

report_lines = []

report_lines.append("="*70)
report_lines.append("BIOMARKER PRIORITIZATION REPORT")
report_lines.append("Day 5: Target/Biomarker Analysis")
report_lines.append("="*70)
report_lines.append("")

report_lines.append("SCORING METHODOLOGY")
report_lines.append("-" * 70)
report_lines.append("Multi-criteria scoring system (0-100 points):")
report_lines.append("  1. Statistical Strength (0-20): Fold change + P-value")
report_lines.append("  2. Mitochondrial Relevance (0-15): Complex assignment")
report_lines.append("  3. Network Centrality (0-15): Hub protein status")
report_lines.append("  5. Disease Association (0-20): Known Leigh/mitochondrial genes")
report_lines.append("  6. Druggability (0-15): Therapeutic potential")
report_lines.append("")

report_lines.append("PRIORITY TIERS")
report_lines.append("-" * 70)
report_lines.append(f"  Tier 1 (≥70 points): {tier_counts.get('Tier 1: High Priority', 0)} candidates - Immediate validation")
report_lines.append(f"  Tier 2 (50-69 points): {tier_counts.get('Tier 2: Medium Priority', 0)} candidates - Secondary targets")
report_lines.append(f"  Tier 3 (30-49 points): {tier_counts.get('Tier 3: Low Priority', 0)} candidates - Exploratory")
report_lines.append(f"  Tier 4 (<30 points): {tier_counts.get('Tier 4: Exploratory', 0)} candidates - Low priority")
report_lines.append("")

report_lines.append("TOP 20 BIOMARKER CANDIDATES")
report_lines.append("-" * 70)
report_lines.append(f"{'Rank':<6}{'Gene':<15}{'Score':<8}{'Tier':<25}{'logFC':<8}{'Mito':<6}{'Known':<7}")
report_lines.append("-" * 70)

for _, row in biomarker_ranking.head(20).iterrows():
    tier_short = row['Priority_Tier'].split(':')[0]
    mito_flag = "✓" if row['Is_Mitochondrial'] else "-"
    known_flag = "✓" if row['Disease_Score'] > 0 else "-"
    report_lines.append(
        f"{int(row['Rank']):<6}{row['Primary_Gene']:<15}{row['Total_Score']:<8.1f}"
        f"{tier_short:<25}{row['logFC']:<8.2f}{mito_flag:<6}{known_flag:<7}"
    )

report_lines.append("")

# Mitochondrial Complex Analysis
report_lines.append("MITOCHONDRIAL COMPLEX DISTRIBUTION")
report_lines.append("-" * 70)

complex_genes = {
    'Complex I': [g for g in biomarker_ranking['Primary_Gene'] if g in complex_I],
    'Complex II': [g for g in biomarker_ranking['Primary_Gene'] if g in complex_II],
    'Complex III': [g for g in biomarker_ranking['Primary_Gene'] if g in complex_III],
    'Complex IV': [g for g in biomarker_ranking['Primary_Gene'] if g in complex_IV],
    'Complex V': [g for g in biomarker_ranking['Primary_Gene'] if g in complex_V]
}

for complex_name, genes in complex_genes.items():
    if genes:
        report_lines.append(f"  {complex_name}: {len(genes)} genes")
        report_lines.append(f"    {', '.join(genes[:10])}")
        if len(genes) > 10:
            report_lines.append(f"    ... and {len(genes) - 10} more")
    else:
        report_lines.append(f"  {complex_name}: No genes found")

report_lines.append("")

# Known Leigh Genes
if len(known_leigh_found) > 0:
    report_lines.append("VALIDATED KNOWN LEIGH SYNDROME GENES")
    report_lines.append("-" * 70)
    for _, row in known_leigh_found.head(10).iterrows():
        report_lines.append(
            f"  {row['Primary_Gene']:15s} Rank: {int(row['Rank']):3d}, "
            f"Score: {row['Total_Score']:.1f}, logFC: {row['logFC']:.2f}"
        )
    report_lines.append("")

# Novel Candidates
novel_candidates = biomarker_ranking[
    (biomarker_ranking['Disease_Score'] == 0) & 
    (biomarker_ranking['Total_Score'] >= 50)
]

if len(novel_candidates) > 0:
    report_lines.append("NOVEL BIOMARKER CANDIDATES (Not Previously Linked to Leigh)")
    report_lines.append("-" * 70)
    for _, row in novel_candidates.head(10).iterrows():
        report_lines.append(
            f"  {row['Primary_Gene']:15s} Score: {row['Total_Score']:.1f}, "
            f"logFC: {row['logFC']:.2f}, Tier: {row['Priority_Tier'].split(':')[0]}"
        )
    report_lines.append("")

# Druggable Targets
druggable_targets = biomarker_ranking[biomarker_ranking['Druggability_Score'] >= 8]

if len(druggable_targets) > 0:
    report_lines.append("POTENTIALLY DRUGGABLE TARGETS")
    report_lines.append("-" * 70)
    for _, row in druggable_targets.head(15).iterrows():
        report_lines.append(
            f"  {row['Primary_Gene']:15s} Druggability: {row['Druggability_Score']:.0f}/15, "
            f"Total Score: {row['Total_Score']:.1f}, Rank: {int(row['Rank'])}"
        )
    report_lines.append("")

# Recommendations
report_lines.append("RECOMMENDATIONS FOR VALIDATION")
report_lines.append("="*70)
report_lines.append("")
report_lines.append("IMMEDIATE PRIORITIES (for experimental validation):")
report_lines.append("")

# Get top tier 1 candidates
tier1_targets = biomarker_ranking[biomarker_ranking['Priority_Tier'] == 'Tier 1: High Priority']

if len(tier1_targets) > 0:
    report_lines.append("1. Western Blot Validation:")
    report_lines.append(f"   Recommended targets: {', '.join(tier1_targets.head(5)['Primary_Gene'].tolist())}")
    report_lines.append("")
    
    report_lines.append("2. qPCR Confirmation:")
    report_lines.append(f"   Top mRNA targets: {', '.join(tier1_targets.head(8)['Primary_Gene'].tolist())}")
    report_lines.append("")
    
    report_lines.append("3. Immunohistochemistry:")
    mito_tier1 = tier1_targets[tier1_targets['Is_Mitochondrial']]
    if len(mito_tier1) > 0:
        report_lines.append(f"   Mitochondrial markers: {', '.join(mito_tier1.head(5)['Primary_Gene'].tolist())}")
    report_lines.append("")
    
    report_lines.append("4. Functional Assays:")
    report_lines.append("   - Respiratory chain complex activity")
    report_lines.append("   - ATP production assays")
    report_lines.append("   - Mitochondrial membrane potential")
    report_lines.append("")

# Statistical Summary
report_lines.append("STATISTICAL SUMMARY")
report_lines.append("-" * 70)
report_lines.append(f"  Total candidates evaluated: {len(biomarker_ranking)}")
report_lines.append(f"  Mitochondrial proteins: {biomarker_ranking['Is_Mitochondrial'].sum()} ({biomarker_ranking['Is_Mitochondrial'].mean()*100:.1f}%)")
report_lines.append(f"  Known disease genes: {(biomarker_ranking['Disease_Score'] > 0).sum()}")
report_lines.append(f"  Novel candidates (Tier 1-2): {len(novel_candidates)}")
report_lines.append(f"  Druggable targets: {len(druggable_targets)}")
report_lines.append(f"  Average score: {biomarker_ranking['Total_Score'].mean():.1f}")
report_lines.append(f"  Median score: {biomarker_ranking['Total_Score'].median():.1f}")
report_lines.append("")

report_lines.append("OUTPUT FILES GENERATED")
report_lines.append("-" * 70)
report_lines.append("  - biomarker_ranking_full.csv           (All candidates)")
report_lines.append("  - biomarker_top50.csv                  (Top 50 targets)")
report_lines.append("  - biomarker_tier1_high_priority.csv    (Immediate targets)")
report_lines.append("  - known_leigh_genes_validated.csv      (Known genes found)")
report_lines.append("  - biomarker_prioritization_summary.png (4-panel figure)")
report_lines.append("  - biomarker_heatmap_top30.png         (Score heatmap)")
report_lines.append("  - biomarker_volcano_prioritized.png    (Prioritized volcano)")
report_lines.append("  - biomarker_report.txt                 (This report)")
report_lines.append("")

report_lines.append("="*70)
report_lines.append("END OF REPORT")
report_lines.append("="*70)

# Save report
report_text = "\n".join(report_lines)
with open(output_dir / 'biomarker_report.txt', 'w') as f:
    f.write(report_text)

print(report_text)
print(f"\n✓ Saved: {output_dir}/biomarker_report.txt")

# ==============================================================================
# 7. CREATE VALIDATION CHECKLIST
# ==============================================================================

print("\n[Step 7] Creating validation checklist for wet lab...\n")

validation_checklist = []

validation_checklist.append("="*70)
validation_checklist.append("EXPERIMENTAL VALIDATION CHECKLIST")
validation_checklist.append("Leigh Syndrome Biomarker Candidates")
validation_checklist.append("="*70)
validation_checklist.append("")

# Phase 1: Tier 1 targets
tier1_list = biomarker_ranking[biomarker_ranking['Priority_Tier'] == 'Tier 1: High Priority']

if len(tier1_list) > 0:
    validation_checklist.append("PHASE 1: HIGH PRIORITY VALIDATION")
    validation_checklist.append("-" * 70)
    validation_checklist.append("")
    
    validation_checklist.append("A. Western Blot Analysis")
    validation_checklist.append("   Target proteins (Top 10):")
    for i, (_, row) in enumerate(tier1_list.head(10).iterrows(), 1):
        validation_checklist.append(f"   □ {i:2d}. {row['Primary_Gene']:15s} - FC: {row['logFC']:+.2f}, P: {row['adj.P.Val']:.2e}")
    validation_checklist.append("")
    validation_checklist.append("   Required antibodies:")
    antibody_list = tier1_list.head(10)['Primary_Gene'].tolist()
    for gene in antibody_list:
        validation_checklist.append(f"   □ Anti-{gene} (recommend Cell Signaling or Abcam)")
    validation_checklist.append("   □ β-Actin or GAPDH (loading control)")
    validation_checklist.append("   □ VDAC1 or COX IV (mitochondrial loading control)")
    validation_checklist.append("")
    
    validation_checklist.append("B. qRT-PCR Validation")
    validation_checklist.append("   Primer design required for:")
    for gene in tier1_list.head(12)['Primary_Gene'].tolist():
        validation_checklist.append(f"   □ {gene}")
    validation_checklist.append("   □ 18S rRNA (housekeeping)")
    validation_checklist.append("   □ ACTB (housekeeping)")
    validation_checklist.append("")
    
    validation_checklist.append("C. Immunofluorescence/IHC")
    mito_targets = tier1_list[tier1_list['Is_Mitochondrial']].head(5)
    validation_checklist.append("   Mitochondrial markers:")
    for _, row in mito_targets.iterrows():
        validation_checklist.append(f"   □ {row['Primary_Gene']} (mitochondrial)")
    validation_checklist.append("   □ MitoTracker (co-staining)")
    validation_checklist.append("   □ DAPI (nuclear staining)")
    validation_checklist.append("")

# Phase 2: Functional assays
validation_checklist.append("PHASE 2: FUNCTIONAL VALIDATION")
validation_checklist.append("-" * 70)
validation_checklist.append("")

validation_checklist.append("A. Mitochondrial Function Assays")
validation_checklist.append("   □ Complex I activity assay")
validation_checklist.append("   □ Complex II activity assay")
validation_checklist.append("   □ Complex IV (COX) activity assay")
validation_checklist.append("   □ Complex V (ATP synthase) activity assay")
validation_checklist.append("   □ Citrate synthase activity (normalization)")
validation_checklist.append("")

validation_checklist.append("B. Cellular Bioenergetics")
validation_checklist.append("   □ ATP production assay (luminescence-based)")
validation_checklist.append("   □ Oxygen consumption rate (Seahorse/Oroboros)")
validation_checklist.append("   □ Mitochondrial membrane potential (TMRE/JC-1)")
validation_checklist.append("   □ ROS production (MitoSOX)")
validation_checklist.append("")

validation_checklist.append("C. Cell/Tissue Models")
validation_checklist.append("   □ Patient-derived fibroblasts (if available)")
validation_checklist.append("   □ Patient muscle biopsies (archived)")
validation_checklist.append("   □ Induced pluripotent stem cells (iPSCs)")
validation_checklist.append("   □ Neuronal differentiation from iPSCs")
validation_checklist.append("")

# Phase 3: Clinical validation
validation_checklist.append("PHASE 3: CLINICAL BIOMARKER VALIDATION")
validation_checklist.append("-" * 70)
validation_checklist.append("")

validation_checklist.append("A. Sample Requirements")
validation_checklist.append("   Patient cohort:")
validation_checklist.append("   □ Leigh syndrome patients (n≥20)")
validation_checklist.append("   □ Age-matched controls (n≥20)")
validation_checklist.append("   □ Other mitochondrial disease controls (n≥10)")
validation_checklist.append("")
validation_checklist.append("   Sample types:")
validation_checklist.append("   □ Plasma/serum (biomarker detection)")
validation_checklist.append("   □ Muscle biopsy (respiratory chain analysis)")
validation_checklist.append("   □ CSF (neurological biomarkers)")
validation_checklist.append("")

validation_checklist.append("B. Biomarker Measurements")
validation_checklist.append("   □ ELISA development for top 5 candidates")
validation_checklist.append("   □ Multi-plex immunoassay (Luminex)")
validation_checklist.append("   □ Mass spectrometry validation")
validation_checklist.append("")

# Sample preparation
validation_checklist.append("SAMPLE PREPARATION PROTOCOLS")
validation_checklist.append("-" * 70)
validation_checklist.append("")

validation_checklist.append("For Western Blot:")
validation_checklist.append("  □ Protein extraction (RIPA or mitochondrial isolation)")
validation_checklist.append("  □ BCA/Bradford protein quantification")
validation_checklist.append("  □ SDS-PAGE (4-12% gradient gels)")
validation_checklist.append("  □ PVDF membrane transfer")
validation_checklist.append("  □ Blocking (5% milk or BSA)")
validation_checklist.append("  □ Primary antibody (1:1000, overnight 4°C)")
validation_checklist.append("  □ Secondary antibody (1:5000, 1h RT)")
validation_checklist.append("  □ ECL detection")
validation_checklist.append("  □ Densitometry analysis (ImageJ)")
validation_checklist.append("")

validation_checklist.append("For qPCR:")
validation_checklist.append("  □ RNA extraction (TRIzol or RNeasy)")
validation_checklist.append("  □ RNA quality check (260/280 ratio ≥1.8)")
validation_checklist.append("  □ DNase treatment")
validation_checklist.append("  □ cDNA synthesis (1-2 μg RNA)")
validation_checklist.append("  □ qPCR (SYBR Green or TaqMan)")
validation_checklist.append("  □ ΔΔCt analysis (normalized to housekeeping)")
validation_checklist.append("")

# Statistical considerations
validation_checklist.append("STATISTICAL CONSIDERATIONS")
validation_checklist.append("-" * 70)
validation_checklist.append("")
validation_checklist.append("Sample Size Calculation:")
validation_checklist.append("  □ Power analysis (α=0.05, β=0.2, effect size from proteomics)")
validation_checklist.append("  □ Minimum n=20 per group recommended")
validation_checklist.append("")
validation_checklist.append("Statistical Tests:")
validation_checklist.append("  □ Normality test (Shapiro-Wilk)")
validation_checklist.append("  □ T-test or Mann-Whitney U (depending on distribution)")
validation_checklist.append("  □ Multiple testing correction (Bonferroni or FDR)")
validation_checklist.append("  □ ROC curve analysis (diagnostic potential)")
validation_checklist.append("")

# Budget estimate
validation_checklist.append("ESTIMATED BUDGET (USD)")
validation_checklist.append("-" * 70)
validation_checklist.append("  Western Blot (10 targets):")
validation_checklist.append("    - Antibodies: $3,000-5,000")
validation_checklist.append("    - Reagents: $500-1,000")
validation_checklist.append("")
validation_checklist.append("  qPCR (12 targets):")
validation_checklist.append("    - Primers: $300-600")
validation_checklist.append("    - Master mix: $500-800")
validation_checklist.append("")
validation_checklist.append("  Functional Assays:")
validation_checklist.append("    - Complex activity kits: $2,000-3,000")
validation_checklist.append("    - ATP assay: $500-800")
validation_checklist.append("    - Seahorse XF analysis: $3,000-5,000")
validation_checklist.append("")
validation_checklist.append("  Total Estimated Cost: $10,000-15,000")
validation_checklist.append("")

validation_checklist.append("="*70)
validation_checklist.append("END OF VALIDATION CHECKLIST")
validation_checklist.append("="*70)

# Save validation checklist
validation_text = "\n".join(validation_checklist)
with open(output_dir / 'validation_checklist.txt', 'w') as f:
    f.write(validation_text)

print("✓ Saved: biomarker_results/validation_checklist.txt\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

print("="*70)
print("DAY 5 BIOMARKER PRIORITIZATION COMPLETE")
print("="*70)
print()
print(f"✓ {len(biomarker_ranking)} candidates evaluated")
print(f"✓ {len(tier1_list)} high-priority targets identified")
print(f"✓ {len(novel_candidates)} novel candidates discovered")
print(f"✓ {len(druggable_targets)} druggable targets found")
print()
print("Output Files Generated:")
print("  1. biomarker_ranking_full.csv")
print("  2. biomarker_top50.csv")
print("  3. biomarker_tier1_high_priority.csv")
print("  4. known_leigh_genes_validated.csv")
print("  5. biomarker_report.txt")
print("  6. validation_checklist.txt")
print("  7. biomarker_prioritization_summary.png")
print("  8. biomarker_heatmap_top30.png")
print("  9. biomarker_volcano_prioritized.png")
print()
print("Next Steps:")
print("  1. Review biomarker_report.txt for detailed findings")
print("  2. Follow validation_checklist.txt for wet-lab experiments")
print("  3. Prioritize Tier 1 candidates for immediate validation")
print("  4. Compare with published Leigh syndrome studies")
print("  5. Prepare manuscript figures from generated plots")
print()
print("="*70)
print("Ready for experimental validation and manuscript preparation!")
print("="*70)