## data-raw/example_PA.R
## Run this script once to regenerate the example datasets.
## Requires: usethis

# =============================================================================
# geneset_list
# A named list of 40 curated gene sets with realistic KEGG / HALLMARK / GO
# naming conventions and real human gene symbols, grouped into four biological
# themes so that geneset_similarity() + do_clust() + get_network_communities()
# produce meaningful clustering results.
#
# Themes:
#   1. Apoptosis & cell death         (8 sets)
#   2. Cell cycle & DNA damage        (8 sets)
#   3. Immune response & inflammation (12 sets)
#   4. Metabolism                     (12 sets)
# =============================================================================

geneset_list <- list(

  # ── 1. Apoptosis & cell death ─────────────────────────────────────────────
  KEGG_APOPTOSIS = c(
    "TP53", "BCL2", "BCL2L1", "BAX", "BAD", "BID", "CASP3", "CASP8",
    "CASP9", "CYCS", "APAF1", "FADD", "FAS", "TNFRSF10A", "TNFRSF10B",
    "MCL1", "XIAP", "DIABLO"
  ),
  HALLMARK_APOPTOSIS = c(
    "TP53", "BCL2", "BAX", "CASP3", "CASP9", "CYCS", "APAF1", "MCL1",
    "BCL2L1", "BID", "PMAIP1", "BBC3", "CASP7", "DFFA", "DFFB"
  ),
  GO_INTRINSIC_APOPTOSIS = c(
    "BAX", "BCL2", "BCL2L1", "BID", "CASP9", "CYCS", "APAF1", "DIABLO",
    "SMAC", "HtrA2", "MCL1", "PMAIP1", "BBC3", "BOK"
  ),
  GO_EXTRINSIC_APOPTOSIS = c(
    "FAS", "FADD", "CASP8", "CASP3", "TNFRSF10A", "TNFRSF10B", "TRADD",
    "RIPK1", "CFLAR", "BID", "CASP10"
  ),
  HALLMARK_P53_PATHWAY = c(
    "TP53", "MDM2", "CDKN1A", "BAX", "PUMA", "NOXA", "GADD45A", "BBC3",
    "PMAIP1", "SIAH1", "BTG2", "SESN1", "SESN2", "TIGAR"
  ),
  KEGG_P53_SIGNALING = c(
    "TP53", "MDM2", "MDM4", "CDKN1A", "GADD45A", "BAX", "BBC3", "PMAIP1",
    "SIAH1", "ATM", "CHEK2", "CDKN2A", "RB1", "CCND1"
  ),
  GO_REGULATION_OF_APOPTOSIS = c(
    "BCL2", "BCL2L1", "MCL1", "XIAP", "BIRC2", "BIRC3", "CFLAR", "TP53",
    "MDM2", "BAX", "BAD", "BID", "CASP3", "CASP8", "DIABLO"
  ),
  GO_MITOCHONDRIAL_OUTER_MEMBRANE_PERMEABILIZATION = c(
    "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "BID", "BIM", "PUMA",
    "NOXA", "VDAC1", "VDAC2", "CYCS"
  ),

  # ── 2. Cell cycle & DNA damage ────────────────────────────────────────────
  KEGG_CELL_CYCLE = c(
    "CDK1", "CDK2", "CDK4", "CDK6", "CCNA1", "CCNA2", "CCNB1", "CCND1",
    "CCNE1", "RB1", "E2F1", "TP53", "CDKN1A", "CDKN2A", "ATM", "CHEK1",
    "CHEK2", "WEE1", "CDC25A", "CDC25C"
  ),
  HALLMARK_E2F_TARGETS = c(
    "E2F1", "E2F2", "E2F3", "CCNA2", "CCNE1", "CDK2", "PCNA", "MCM2",
    "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "RRM1", "RRM2", "TYMS"
  ),
  HALLMARK_G2M_CHECKPOINT = c(
    "CDK1", "CCNB1", "CCNB2", "CDC20", "BUB1", "BUB1B", "MAD2L1",
    "AURKA", "AURKB", "PLK1", "WEE1", "CHEK1", "CDC25C", "BRCA1"
  ),
  GO_DNA_DAMAGE_RESPONSE = c(
    "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2", "RAD51", "H2AFX",
    "MDC1", "RNF8", "RNF168", "PALB2", "FANCD2", "TP53", "CDKN1A"
  ),
  KEGG_DNA_REPLICATION = c(
    "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "RFC1",
    "RFC2", "RFC3", "RFC4", "RFC5", "POLA1", "POLD1", "POLE", "LIG1"
  ),
  GO_MITOTIC_CELL_CYCLE = c(
    "CDK1", "CCNB1", "CDC20", "APC", "MAD2L1", "BUB1", "BUB1B",
    "AURKA", "AURKB", "PLK1", "ESPL1", "SECURIN", "SMC1A", "SMC3"
  ),
  HALLMARK_MYC_TARGETS_V1 = c(
    "MYC", "CDK4", "CCND1", "E2F1", "PCNA", "MCM2", "NPM1", "RCC1",
    "POLD1", "RFC2", "TERT", "LDHA", "ENO1", "PKM", "TKT"
  ),
  KEGG_HOMOLOGOUS_RECOMBINATION = c(
    "BRCA1", "BRCA2", "RAD51", "RAD52", "PALB2", "RBBP8", "MRE11",
    "RAD50", "NBN", "ATM", "FANCD2", "XRCC2", "XRCC3"
  ),

  # ── 3. Immune response & inflammation ─────────────────────────────────────
  HALLMARK_INFLAMMATORY_RESPONSE = c(
    "IL1B", "IL6", "TNF", "CXCL8", "CCL2", "PTGS2", "NF-kB1", "RELA",
    "NFKBIA", "ICAM1", "VCAM1", "SELE", "IL1A", "CXCL1", "CXCL2"
  ),
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = c(
    "TNF", "TNFRSF1A", "TRADD", "TRAF2", "RIPK1", "NFKB1", "RELA",
    "NFKBIA", "NFKBIB", "IKBKA", "IKBKB", "IKBKG", "IL6", "IL8", "ICAM1"
  ),
  KEGG_NF_KAPPA_B_SIGNALING = c(
    "NFKB1", "NFKB2", "RELA", "RELB", "REL", "NFKBIA", "NFKBIB",
    "IKBKA", "IKBKB", "IKBKG", "TNF", "IL1B", "LTA", "CD40LG"
  ),
  GO_CYTOKINE_MEDIATED_SIGNALING = c(
    "IL6", "IL6R", "JAK1", "JAK2", "STAT1", "STAT3", "IL2", "IL4",
    "IL10", "IL12A", "IFNG", "IFNA1", "IRF1", "IRF3", "SOCS1", "SOCS3"
  ),
  HALLMARK_INTERFERON_GAMMA_RESPONSE = c(
    "IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1",
    "CIITA", "HLA-DRA", "HLA-A", "TAP1", "TAP2", "B2M", "PSMB9", "CXCL10"
  ),
  HALLMARK_INTERFERON_ALPHA_RESPONSE = c(
    "IFNA1", "IFNAR1", "IFNAR2", "JAK1", "TYK2", "STAT1", "STAT2",
    "IRF3", "IRF7", "ISG15", "MX1", "OAS1", "OAS2", "IFIT1", "IFIT3"
  ),
  KEGG_JAK_STAT_SIGNALING = c(
    "JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT2", "STAT3", "STAT5A",
    "STAT5B", "IL2", "IL6", "IFNG", "IFNA1", "SOCS1", "SOCS3", "PTPN11"
  ),
  GO_T_CELL_ACTIVATION = c(
    "CD3E", "CD3G", "CD247", "ZAP70", "LCK", "FYN", "LAT", "PLCG1",
    "NFATC1", "NFATC2", "IL2", "CD28", "CD80", "CD86", "CTLA4", "PDCD1"
  ),
  KEGG_TOLL_LIKE_RECEPTOR_SIGNALING = c(
    "TLR1", "TLR2", "TLR4", "TLR9", "MYD88", "TRIF", "IRAK1", "IRAK4",
    "TRAF6", "NFKB1", "IRF3", "IRF7", "TNF", "IL6", "IL12A"
  ),
  GO_COMPLEMENT_ACTIVATION = c(
    "C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B", "C5", "CFB",
    "CFD", "CFH", "CFI", "CR1", "CD55", "CD59"
  ),
  HALLMARK_IL6_JAK_STAT3_SIGNALING = c(
    "IL6", "IL6R", "IL6ST", "JAK1", "JAK2", "STAT3", "SOCS1", "SOCS3",
    "IL10", "IL11", "IL22", "IL31RA", "OSMR", "LIFR"
  ),
  KEGG_CHEMOKINE_SIGNALING = c(
    "CCL2", "CCL5", "CXCL8", "CXCL10", "CCR2", "CCR5", "CXCR3", "CXCR4",
    "GNB1", "GNG2", "PLCB1", "PIK3CA", "AKT1", "MAPK1", "STAT3"
  ),

  # ── 4. Metabolism ─────────────────────────────────────────────────────────
  HALLMARK_GLYCOLYSIS = c(
    "HK1", "HK2", "GPI", "PFKL", "PFKM", "ALDOA", "TPI1", "GAPDH",
    "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "SLC2A1", "SLC2A3"
  ),
  KEGG_GLYCOLYSIS_GLUCONEOGENESIS = c(
    "HK1", "HK2", "GPI", "PFKL", "ALDOA", "TPI1", "GAPDH", "PGK1",
    "ENO1", "PKM", "LDHA", "PCK1", "PCK2", "G6PC", "FBP1", "FBP2"
  ),
  GO_GLUCONEOGENESIS = c(
    "PCK1", "PCK2", "FBP1", "FBP2", "G6PC", "ALDOB", "ENO1", "GAPDH",
    "PGK1", "PGAM1", "MDH1", "MDH2", "GOT1", "GOT2"
  ),
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = c(
    "NDUFA1", "NDUFA2", "NDUFB1", "SDHA", "SDHB", "UQCRC1", "UQCRC2",
    "COX4I1", "COX5A", "ATP5A1", "ATP5B", "ATP5C1", "CYCS", "SCO1"
  ),
  KEGG_CITRATE_CYCLE_TCA = c(
    "CS", "ACO1", "ACO2", "IDH1", "IDH2", "OGDH", "SUCLA2", "SUCLG1",
    "SDHA", "SDHB", "FH", "MDH1", "MDH2", "DLST", "DLD"
  ),
  GO_FATTY_ACID_BETA_OXIDATION = c(
    "HADHA", "HADHB", "ACADL", "ACADM", "ACADS", "ACADVL", "ECHS1",
    "HADH", "CPT1A", "CPT1B", "CPT2", "SLC25A20", "ACSL1", "ACSL3"
  ),
  KEGG_FATTY_ACID_METABOLISM = c(
    "FASN", "ACACA", "ACACB", "ACSL1", "ACSL3", "ACSL4", "HADHA",
    "HADHB", "ACADL", "ACADM", "CPT1A", "CPT2", "ELOVL1", "ELOVL6"
  ),
  HALLMARK_FATTY_ACID_METABOLISM = c(
    "FASN", "ACACA", "ACSL1", "HADHA", "ACADL", "CPT1A", "PPARA",
    "PPARG", "RXRA", "FABP1", "FABP4", "ACOX1", "EHHADH", "HSD17B4"
  ),
  KEGG_PENTOSE_PHOSPHATE_PATHWAY = c(
    "G6PD", "PGD", "RPIA", "TALDO1", "TKT", "RPE", "RBKS", "PGLS",
    "H6PD", "DERA", "PRPS1", "PRPS2"
  ),
  KEGG_PURINE_METABOLISM = c(
    "ADSS", "ADSL", "ATIC", "PFAS", "GART", "PPAT", "PAICS", "ADA",
    "AMPD1", "DGUOK", "HPRT1", "IMPDH1", "IMPDH2", "GMPS"
  ),
  KEGG_AMINO_SUGAR_NUCLEOTIDE_SUGAR_METABOLISM = c(
    "GFPT1", "GFPT2", "UAP1", "GNPNAT1", "PGM3", "NAGK", "CMAS",
    "SLC35A1", "UGP2", "UGDH", "B4GALT1", "ST6GAL1"
  ),
  HALLMARK_MTORC1_SIGNALING = c(
    "MTOR", "RPTOR", "RPS6KB1", "RPS6KB2", "EIF4EBP1", "EIF4E", "AKT1",
    "TSC1", "TSC2", "RHEB", "DEPTOR", "MLST8", "PRAS40", "HIF1A", "MYC"
  )
)

# =============================================================================
# camera_results
# A data.frame simulating CAMERA output for a differential expression analysis.
# Columns: GeneSet, Direction, PValue, FDR.
# ~60 % of gene sets are significant (FDR < 0.05) so clustering has enough
# sets to work with.
# =============================================================================

set.seed(174)

n_sets   <- length(geneset_list)
set_names <- names(geneset_list)

# Simulate realistic p-value distribution
raw_p <- c(
  # Apoptosis theme — mostly significant
  runif(8,  min = 1e-6, max = 0.02),
  # Cell cycle theme — mostly significant
  runif(8,  min = 1e-5, max = 0.03),
  # Immune theme — mix
  c(runif(7, min = 1e-4, max = 0.04), runif(5, min = 0.06, max = 0.50)),
  # Metabolism theme — mix
  c(runif(7, min = 1e-3, max = 0.04), runif(5, min = 0.10, max = 0.80))
)

fdr_vals  <- p.adjust(raw_p, method = "BH")
direction <- ifelse(
  runif(n_sets) > 0.4, "Up", "Down"
)

camera_results <- data.frame(
  GeneSet   = set_names,
  Direction = direction,
  PValue    = raw_p,
  FDR       = fdr_vals,
  stringsAsFactors = FALSE
)

# =============================================================================
# Save to data/
# =============================================================================

usethis::use_data(geneset_list,   overwrite = TRUE)
usethis::use_data(camera_results, overwrite = TRUE)

message("Done. Objects saved:")
message("  data/geneset_list.rda   — ", length(geneset_list), " gene sets")
message("  data/camera_results.rda — ", nrow(camera_results), " rows, ",
        sum(camera_results$FDR < 0.05), " sets with FDR < 0.05")
