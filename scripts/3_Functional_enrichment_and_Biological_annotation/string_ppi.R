#!/usr/bin/env Rscript
# ==============================================================================
# STRING Protein-Protein Interaction Network Analysis (Reproducible)
# Using STRINGdb Bioconductor package
# ==============================================================================

# Suppress warnings
options(warn = -1)

cat("======================================================================\n")
cat("STRING PROTEIN-PROTEIN INTERACTION ANALYSIS\n")
cat("======================================================================\n\n")

# ==============================================================================
# 1. LOAD REQUIRED PACKAGES
# ==============================================================================

cat("[Step 1] Loading required packages...\n")

# Check and install required packages
required_packages <- c("STRINGdb", "igraph", "ggplot2", "dplyr", "readr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    if (pkg == "STRINGdb") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install("STRINGdb", update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

cat("✓ All packages loaded\n\n")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

cat("[Step 2] Configuration...\n")

# Input/Output paths
INPUT_FILE <- "differential_results/sig_for_enrichment.tsv"
OUTPUT_DIR <- "new_1_enrichment_results"

# STRING parameters
SPECIES <- 9606  # 9606 = Homo sapiens, 10090 = Mus musculus
SCORE_THRESHOLD <- 400  # Medium confidence (0-1000 scale)
STRING_VERSION <- "12.0"
MAX_GENES <- 200  # Maximum genes to analyze (for performance)

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("  Input file: %s\n", INPUT_FILE))
cat(sprintf("  Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("  Species: %d (Human)\n", SPECIES))
cat(sprintf("  Score threshold: %d (medium confidence)\n", SCORE_THRESHOLD))
cat(sprintf("  STRING version: %s\n\n", STRING_VERSION))

# ==============================================================================
# 3. LOAD DIFFERENTIAL EXPRESSION RESULTS
# ==============================================================================

cat("[Step 3] Loading differential expression results...\n")

# Read significant proteins
sig_results <- read_delim(INPUT_FILE, delim = "\t", show_col_types = FALSE)

cat(sprintf("✓ Loaded %d significant genes\n", nrow(sig_results)))
cat(sprintf("  Up-regulated: %d\n", sum(sig_results$Regulation == "Up")))
cat(sprintf("  Down-regulated: %d\n\n", sum(sig_results$Regulation == "Down")))

# Prepare gene list - ensure proper data frame structure
genes_df <- sig_results %>%
  select(Gene_names, logFC, adj.P.Val, Regulation) %>%
  filter(!is.na(Gene_names), Gene_names != "", Gene_names != "NA")

# Convert tibble to data.frame (STRINGdb requires base R data.frame)
genes_df <- as.data.frame(genes_df)

# Ensure Gene_names is character, not factor
genes_df$Gene_names <- as.character(genes_df$Gene_names)

# Remove any duplicate genes (take first occurrence)
genes_df <- genes_df[!duplicated(genes_df$Gene_names), ]

# Limit to top genes by significance (if needed)
if (nrow(genes_df) > MAX_GENES) {
  cat(sprintf("\n  Limiting to top %d genes by p-value...\n", MAX_GENES))
  genes_df <- genes_df[order(genes_df$adj.P.Val), ]
  genes_df <- genes_df[1:MAX_GENES, ]
}

cat(sprintf("  Final gene list: %d genes for STRING analysis\n", nrow(genes_df)))

# Save the gene list for reproducibility (exactly what STRING will see)
gene_list_for_string <- data.frame(
  Gene = genes_df$Gene_names,
  stringsAsFactors = FALSE
)

write.table(gene_list_for_string$Gene,
            file.path(OUTPUT_DIR, "string_input_genes.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

cat(sprintf("✓ Saved gene list: %s/string_input_genes.txt\n", OUTPUT_DIR))
cat("  (Use this exact list to reproduce in STRING web interface)\n\n")

# ==============================================================================
# 4. INITIALIZE STRING DATABASE
# ==============================================================================

cat("[Step 4] Initializing STRING database...\n")

string_db <- STRINGdb$new(
  version = STRING_VERSION,
  species = SPECIES,
  score_threshold = SCORE_THRESHOLD,
  input_directory = ""
)

cat("✓ STRING database initialized\n\n")

# ==============================================================================
# 5. MAP GENES TO STRING IDs
# ==============================================================================

cat("[Step 5] Mapping genes to STRING identifiers...\n")

# CRITICAL: STRINGdb::map() is very picky about input format
# Ensure proper data.frame structure
genes_for_mapping <- data.frame(
  Gene_names = as.character(genes_df$Gene_names),
  logFC = as.numeric(genes_df$logFC),
  adj.P.Val = as.numeric(genes_df$adj.P.Val),
  Regulation = as.character(genes_df$Regulation),
  stringsAsFactors = FALSE
)

# Remove any rows with NA in Gene_names
genes_for_mapping <- genes_for_mapping[!is.na(genes_for_mapping$Gene_names), ]

cat(sprintf("  Input: %d genes to map\n", nrow(genes_for_mapping)))
cat(sprintf("  Column name: '%s'\n", "Gene_names"))

# Map gene names to STRING IDs
mapped <- tryCatch({
  string_db$map(
    genes_for_mapping,
    "Gene_names",
    removeUnmappedRows = TRUE,
    takeFirst = TRUE
  )
}, error = function(e) {
  cat("\n⚠ Mapping error occurred. Trying alternative approach...\n")
  cat(sprintf("  Error message: %s\n\n", e$message))
  
  # Alternative: Simple data frame with just gene names
  simple_df <- data.frame(
    Gene_names = genes_for_mapping$Gene_names,
    stringsAsFactors = FALSE
  )
  
  mapped_alt <- string_db$map(
    simple_df,
    "Gene_names",
    removeUnmappedRows = FALSE,
    takeFirst = TRUE
  )
  
  # Add back the metadata
  mapped_alt <- merge(
    mapped_alt,
    genes_for_mapping,
    by = "Gene_names",
    all.x = TRUE
  )
  
  # Remove unmapped
  mapped_alt <- mapped_alt[!is.na(mapped_alt$STRING_id), ]
  
  return(mapped_alt)
})

n_mapped <- nrow(mapped)
n_unmapped <- nrow(genes_df) - n_mapped
mapping_rate <- n_mapped / nrow(genes_df) * 100

cat(sprintf("✓ Mapping complete:\n"))
cat(sprintf("  Mapped: %d genes (%.1f%%)\n", n_mapped, mapping_rate))
cat(sprintf("  Unmapped: %d genes\n\n", n_unmapped))

# Check unmapped genes
if (n_unmapped > 0) {
  unmapped_genes <- setdiff(genes_df$Gene_names, mapped$Gene_names)
  cat("  Top unmapped genes:\n")
  cat(paste0("    - ", head(unmapped_genes, 10), collapse = "\n"), "\n\n")
}

# Save mapping results
write.csv(mapped, file.path(OUTPUT_DIR, "string_gene_mapping.csv"), row.names = FALSE)
cat(sprintf("✓ Saved: %s/string_gene_mapping.csv\n\n", OUTPUT_DIR))

# ==============================================================================
# 6. GET PROTEIN-PROTEIN INTERACTIONS
# ==============================================================================

cat("[Step 6] Retrieving protein-protein interactions...\n")

# Get interactions from STRING
interactions <- string_db$get_interactions(mapped$STRING_id)
cat(" Downloading interaction data...\n")

n_interactions <- nrow(interactions)
n_proteins <- length(unique(c(interactions$from, interactions$to)))

cat(sprintf("✓ Retrieved interactions:\n"))
cat(sprintf("  Total interactions: %d\n", n_interactions))
cat(sprintf("  Unique proteins: %d\n", n_proteins))

if (n_interactions > 0) {
  cat(sprintf("  Average combined score: %.3f\n\n", mean(interactions$combined_score)))
  
  # Add gene names to interactions
  interactions_annotated <- interactions %>%
    left_join(mapped %>% select(STRING_id, Gene_names, logFC, Regulation),
              by = c("from" = "STRING_id")) %>%
    rename(gene1 = Gene_names, logFC1 = logFC, regulation1 = Regulation) %>%
    left_join(mapped %>% select(STRING_id, Gene_names, logFC, Regulation),
              by = c("to" = "STRING_id")) %>%
    rename(gene2 = Gene_names, logFC2 = logFC, regulation2 = Regulation)
  
  # Save interactions
  write.csv(interactions_annotated, 
            file.path(OUTPUT_DIR, "string_interactions.csv"), 
            row.names = FALSE)
  cat(sprintf("✓ Saved: %s/string_interactions.csv\n\n", OUTPUT_DIR))
  
} else {
  cat("⚠ No interactions found\n\n")
  quit(save = "no", status = 0)
}

# ==============================================================================
# 7. NETWORK STATISTICS
# ==============================================================================

cat("[Step 7] Calculating network statistics...\n")

# Create igraph network
g <- graph_from_data_frame(
  interactions_annotated[, c("gene1", "gene2", "combined_score")],
  directed = FALSE
)

# Calculate network metrics
degree_vals <- degree(g)
betweenness_vals <- betweenness(g)
closeness_vals <- closeness(g)

# Hub proteins (top 10 by degree)
hub_proteins <- sort(degree_vals, decreasing = TRUE)[1:10]

cat("\nNetwork Statistics:\n")
cat(sprintf("  Nodes: %d\n", vcount(g)))
cat(sprintf("  Edges: %d\n", ecount(g)))
cat(sprintf("  Average degree: %.2f\n", mean(degree_vals)))
cat(sprintf("  Network density: %.4f\n", edge_density(g)))
cat(sprintf("  Connected components: %d\n", components(g)$no))
cat(sprintf("  Clustering coefficient: %.3f\n", transitivity(g)))

cat("\n  Top 10 Hub Proteins (by degree):\n")
for (i in 1:length(hub_proteins)) {
  cat(sprintf("    %2d. %s: %d connections\n", i, names(hub_proteins)[i], hub_proteins[i]))
}
cat("\n")

# Save network statistics
network_stats <- data.frame(
  Gene = V(g)$name,
  Degree = degree_vals,
  Betweenness = betweenness_vals,
  Closeness = closeness_vals,
  stringsAsFactors = FALSE
) %>%
  left_join(mapped %>% select(Gene_names, logFC, adj.P.Val, Regulation),
            by = c("Gene" = "Gene_names")) %>%
  arrange(desc(Degree))

write.csv(network_stats, 
          file.path(OUTPUT_DIR, "string_network_statistics.csv"), 
          row.names = FALSE)
cat(sprintf("✓ Saved: %s/string_network_statistics.csv\n\n", OUTPUT_DIR))

# ==============================================================================
# 8. COMMUNITY DETECTION
# ==============================================================================

cat("[Step 8] Detecting functional modules...\n")

# Detect communities using Louvain algorithm
communities <- cluster_louvain(g)

cat(sprintf("✓ Found %d functional modules\n", length(communities)))
cat(sprintf("  Modularity: %.3f\n\n", modularity(communities)))

# Module membership
module_membership <- data.frame(
  Gene = V(g)$name,
  Module = membership(communities),
  stringsAsFactors = FALSE
) %>%
  left_join(mapped %>% select(Gene_names, logFC, Regulation),
            by = c("Gene" = "Gene_names"))

# Count genes per module
module_sizes <- table(module_membership$Module)
cat("  Module sizes:\n")
for (i in 1:min(10, length(module_sizes))) {
  cat(sprintf("    Module %d: %d genes\n", i, module_sizes[i]))
}
cat("\n")

write.csv(module_membership, 
          file.path(OUTPUT_DIR, "string_modules.csv"), 
          row.names = FALSE)
cat(sprintf("✓ Saved: %s/string_modules.csv\n\n", OUTPUT_DIR))

# ==============================================================================
# 9. FUNCTIONAL ENRICHMENT OF NETWORK
# ==============================================================================

cat("[Step 9] Performing functional enrichment of network...\n")

# Try multiple enrichment categories
enrichment_categories <- c("Process", "Component", "Function", "KEGG", "Pfam")

for (category in enrichment_categories) {
  cat(sprintf("  Trying %s enrichment...\n", category))
  
  # Get enrichment for the network
  enrichment <- string_db$get_enrichment(mapped$STRING_id, category = category)
  
  if (!is.null(enrichment) && nrow(enrichment) > 0) {
    filename <- sprintf("string_enrichment_%s.csv", tolower(category))
    write.csv(enrichment, file.path(OUTPUT_DIR, filename), row.names = FALSE)
    cat(sprintf("✓ Found %d terms, saved to %s\n", nrow(enrichment), filename))
    
    # Show top 5 for first successful category
    if (category == enrichment_categories[1]) {
      cat("\n  Top 5 enriched terms:\n")
      top_terms <- head(enrichment[order(enrichment$p_value), ], 5)
      for (i in seq_len(nrow(top_terms))) {
        cat(sprintf("    %d. %s (p=%.2e)\n", 
                    i, 
                    substr(top_terms$term_description[i], 1, 50),
                    top_terms$p_value[i]))
      }
    }
  } else {
    cat(sprintf("    ⚠ No %s enrichment found\n", category))
  }
}
cat("\n")

# ==============================================================================
# 10. VISUALIZATION
# ==============================================================================

cat("[Step 10] Creating network visualizations...\n")

# Set seed for reproducibility
set.seed(42)

# --- Plot 1: Full Network with Communities ---
png(file.path(OUTPUT_DIR, "string_network_communities.png"), 
    width = 3000, height = 3000, res = 300)

# Color by module
module_colors <- rainbow(length(communities))
vertex_colors <- module_colors[membership(communities)]

# Size by degree
vertex_sizes <- scales::rescale(degree(g), to = c(5, 20))

layout_fr <- layout_with_fr(g)

plot(communities, g,
     layout = layout_fr,
     vertex.size = vertex_sizes,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     edge.width = 0.5,
     edge.color = rgb(0.5, 0.5, 0.5, 0.3),
     main = "STRING Protein-Protein Interaction Network\n(Colored by Functional Module)")

dev.off()
cat("✓ Saved: new_1_enrichment_results/string_network_communities.png\n")

# --- Plot 2: Network colored by Regulation ---
png(file.path(OUTPUT_DIR, "string_network_regulation.png"), 
    width = 3000, height = 3000, res = 300)

# Color by regulation
reg_colors <- ifelse(module_membership$Regulation == "Up", 
                     "#d62728",  # Red for up
                     ifelse(module_membership$Regulation == "Down",
                            "#1f77b4",  # Blue for down
                            "gray"))    # Gray for not significant

plot(g,
     layout = layout_fr,
     vertex.color = reg_colors,
     vertex.size = vertex_sizes,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     edge.width = 0.5,
     edge.color = rgb(0.5, 0.5, 0.5, 0.3),
     main = "STRING Network: Protein Regulation\n(Red=Up, Blue=Down)")

legend("topright",
       legend = c("Up-regulated", "Down-regulated"),
       col = c("#d62728", "#1f77b4"),
       pch = 19,
       pt.cex = 2,
       bty = "n")

dev.off()
cat("✓ Saved: new_1_enrichment_results/string_network_regulation.png\n")

# --- Plot 3: Hub Protein Subnetwork (Top 30) ---
if (n_proteins > 30) {
  top_hubs <- names(sort(degree_vals, decreasing = TRUE)[1:30])
  g_sub <- induced_subgraph(g, top_hubs)
  
  png(file.path(OUTPUT_DIR, "string_network_hubs.png"), 
      width = 3000, height = 3000, res = 300)
  
  reg_colors_sub <- reg_colors[V(g)$name %in% top_hubs]
  vertex_sizes_sub <- scales::rescale(degree(g_sub), to = c(8, 25))
  
  plot(g_sub,
       layout = layout_with_fr(g_sub),
       vertex.color = reg_colors_sub,
       vertex.size = vertex_sizes_sub,
       vertex.label.cex = 1.0,
       vertex.label.color = "black",
       vertex.label.family = "sans",
       edge.width = 1,
       edge.color = rgb(0.3, 0.3, 0.3, 0.5),
       main = "Top 30 Hub Proteins\n(Size = Connectivity)")
  
  dev.off()
  cat("✓ Saved: new_1_enrichment_results/string_network_hubs.png\n")
}

# --- Plot 4: Network Statistics Summary ---
png(file.path(OUTPUT_DIR, "string_network_statistics_plot.png"), 
    width = 3600, height = 2400, res = 300)

par(mfrow = c(2, 3), mar = c(5, 4, 4, 2))

# Degree distribution
hist(degree_vals, breaks = 30, 
     col = "steelblue", 
     main = "Degree Distribution",
     xlab = "Degree", 
     ylab = "Frequency")

# Betweenness distribution
hist(log10(betweenness_vals + 1), breaks = 30,
     col = "coral",
     main = "Betweenness Distribution",
     xlab = "log10(Betweenness + 1)",
     ylab = "Frequency")

# Degree vs Betweenness
plot(degree_vals, log10(betweenness_vals + 1),
     pch = 19, col = rgb(0, 0, 0, 0.3),
     xlab = "Degree",
     ylab = "log10(Betweenness + 1)",
     main = "Degree vs Betweenness")

# Module size distribution
barplot(sort(table(module_membership$Module), decreasing = TRUE),
        col = "purple",
        main = "Module Sizes",
        xlab = "Module",
        ylab = "Number of Proteins")

# Edge weight distribution
hist(interactions$combined_score, breaks = 30,
     col = "darkgreen",
     main = "Interaction Score Distribution",
     xlab = "Combined Score",
     ylab = "Frequency")

dev.off()
cat("✓ Saved: new_1_enrichment_results/string_network_statistics_plot.png\n\n")

# ==============================================================================
# 11. EXPORT FOR CYTOSCAPE
# ==============================================================================

cat("[Step 11] Exporting network for Cytoscape...\n")

# Node attributes
node_attributes <- network_stats %>%
  select(Gene, Degree, Betweenness, Closeness, logFC, adj.P.Val, Regulation)

write.table(node_attributes,
            file.path(OUTPUT_DIR, "cytoscape_nodes.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Edge list
edge_list <- interactions_annotated %>%
  select(gene1, gene2, combined_score)

write.table(edge_list,
            file.path(OUTPUT_DIR, "cytoscape_edges.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("✓ Saved Cytoscape files:\n")
cat("  - cytoscape_nodes.txt\n")
cat("  - cytoscape_edges.txt\n\n")

cat("  Import to Cytoscape:\n")
cat("    1. File → Import → Network from File\n")
cat("    2. Select cytoscape_edges.txt\n")
cat("    3. File → Import → Table from File\n")
cat("    4. Select cytoscape_nodes.txt\n\n")

# ==============================================================================
# 12. SUMMARY REPORT
# ==============================================================================

cat("======================================================================\n")
cat("STRING ANALYSIS COMPLETE - SUMMARY REPORT\n")
cat("======================================================================\n\n")

cat(sprintf("Input:\n"))
cat(sprintf("  Genes analyzed: %d\n", nrow(genes_df)))
cat(sprintf("  Genes mapped to STRING: %d (%.1f%%)\n", n_mapped, mapping_rate))

cat(sprintf("\nNetwork:\n"))
cat(sprintf("  Proteins: %d\n", n_proteins))
cat(sprintf("  Interactions: %d\n", n_interactions))
cat(sprintf("  Average connections per protein: %.2f\n", mean(degree_vals)))
cat(sprintf("  Network density: %.4f\n", edge_density(g)))

cat(sprintf("\nModules:\n"))
cat(sprintf("  Functional modules detected: %d\n", length(communities)))
cat(sprintf("  Modularity: %.3f\n", modularity(communities)))

cat(sprintf("\nTop Hub Proteins:\n"))
for (i in 1:min(5, length(hub_proteins))) {
  cat(sprintf("  %d. %s (%d connections)\n", i, names(hub_proteins)[i], hub_proteins[i]))
}

cat(sprintf("\nOutput Files:\n"))
cat("  - string_gene_mapping.csv\n")
cat("  - string_interactions.csv\n")
cat("  - string_network_statistics.csv\n")
cat("  - string_modules.csv\n")
cat("  - string_enrichment_*.csv\n")
cat("  - string_network_communities.png\n")
cat("  - string_network_regulation.png\n")
cat("  - string_network_hubs.png\n")
cat("  - string_network_statistics_plot.png\n")
cat("  - cytoscape_nodes.txt\n")
cat("  - cytoscape_edges.txt\n")

cat("\n======================================================================\n")
cat("Analysis complete! Results saved to: new_1_enrichment_results/\n")