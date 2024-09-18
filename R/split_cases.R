########################
# Function split_cases #
########################

#' Obtain 10 exclusive cases from 3 comparisons.
#'
#' When performing differential expression analysis of a study with 3 phenotypes,
#' including the baseline, there are 10 mutually exclusive cases where genes can
#' fall into. This function allows us to obtain these 10 cases and saves them
#' into a list.
#' 
#' @param df.BvsA Data frame comparing the first condition to the baseline.
#' @param df.CvsA Data frame comparing the second condition to the baseline.
#' @param df.BvsC Data frame comparing the two conditions of the study.
#' @param unique_id Column name with the unique identifiers in the data tables (i.e. ensembl for DESeq2, set NAME for GSEA).
#' @param significance_var Variable of significance to filter by (padj for DESeq2 or FDR for GSEA).
#' @param significance_cutoff Cut-off of the significance variable.
#' @param change_var Variable that indicates the direction of the change (i.e. log2FoldChange in DESeq2, NES in GSEA).
#' @param change_cutoff The values of the change variable will be filtered by |change_var| > change_cutoff. Default: 0.
#' @export


split_cases <- function (df.BvsA = NULL, df.CvsA = NULL, df.BvsC = NULL, unique_id = "ensembl",
                         significance_var = "padj", significance_cutoff = 0.25,
                         change_var = "log2FoldChange", change_cutoff = 0)
  
{
  # Set row names as unique identifiers
  if (!all(rownames(df.BvsA) == df.BvsA[, unique_id])) {
    rownames(df.BvsA) <- df.BvsA[, unique_id]
  }
  
  if (!all(rownames(df.CvsA) == df.CvsA[, unique_id])) {
    rownames(df.CvsA) <- df.CvsA[, unique_id]
  }
  
  if (!all(rownames(df.BvsC) == df.BvsC[, unique_id])) {
    rownames(df.BvsC) <- df.BvsC[, unique_id]
  }

  
  # Create subsets by significance
  df.BvsA.sig <- df.BvsA[df.BvsA[, significance_var] < significance_cutoff & !is.na(df.BvsA[, significance_var]), ]
  df.CvsA.sig <- df.CvsA[df.CvsA[, significance_var] < significance_cutoff & !is.na(df.CvsA[, significance_var]), ]
  df.BvsC.sig <- df.BvsC[df.BvsC[, significance_var] < significance_cutoff & !is.na(df.BvsC[, significance_var]), ]
  
  df.BvsA.nsig <- df.BvsA[df.BvsA[, significance_var] >= significance_cutoff | is.na(df.BvsA[, significance_var]), ]
  df.CvsA.nsig <- df.CvsA[df.CvsA[, significance_var] >= significance_cutoff | is.na(df.CvsA[, significance_var]), ]
  df.BvsC.nsig <- df.BvsC[df.BvsC[, significance_var] >= significance_cutoff | is.na(df.BvsC[, significance_var]), ]
  
  
  # Defining cases
  
  ####################### Significant in all comparisons #######################
  
  #############################
  # Case 1: Ladders up/down
  #############################
  # These genes would show a progressively increasing or decreasing expression,
  # particularly useful when comparing conditions over time, intensity or evolution
  
  
  case1.up.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] > change_cutoff, unique_id]))
  
  case1.up <- as.data.frame(df.BvsA.sig[case1.up.genes, ])
  case1.up$trend <- rep("up", nrow(case1.up))
  
  case1.dn.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] < -change_cutoff, unique_id]))
  
  case1.dn <- as.data.frame(df.BvsA.sig[case1.dn.genes, ])
  case1.dn$trend <- rep("dn", nrow(case1.dn))
  
  ##################################
  # Case 2: Stronger in condition 1 
  ##################################
  # Genes with stronger dysregulation in the first condition, which means that
  # they change their expression trend while still being significant
  
  case2.up.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] < -change_cutoff, unique_id]))
  
  case2.up <- as.data.frame(df.BvsA.sig[case2.up.genes, ])
  case2.up$trend <- rep("up", nrow(case2.up))
  
  case2.dn.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] > change_cutoff, unique_id]))
  
  case2.dn <- as.data.frame(df.BvsA.sig[case2.dn.genes, ])
  case2.dn$trend <- rep("dn", nrow(case2.dn))
  
  ##################################
  # Case 3: Stronger in condition 2
  ##################################
  # Genes with stronger dysregulation in the second condition, which means that
  # they change their expression trend while still being significant
  
  case3.up.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] > change_cutoff, unique_id]))
  
  case3.up <- as.data.frame(df.BvsA.sig[case3.up.genes, ])
  case3.up$trend <- rep("up", nrow(case3.up))
  
  case3.dn.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] < -change_cutoff, unique_id]))
  
  case3.dn <- as.data.frame(df.BvsA.sig[case3.dn.genes, ])
  case3.dn$trend <- rep("dn", nrow(case3.dn))
  
  
  ##################### Significant in only two comparisons #####################
  
  ###########################################
  # Case 4: Significant in data 2 and data 3
  ###########################################
  # These genes would be consider markers of the second condition since they are
  # dysregulated compared to the baseline and to the first condition only
  
  case4.up.genes <- Reduce(intersect, list(df.CvsA.sig[df.CvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] > change_cutoff, unique_id],
                                           df.BvsA.nsig[, unique_id]))
  
  case4.up <- as.data.frame(df.CvsA.sig[case4.up.genes, ])
  case4.up$trend <- rep("up", nrow(case4.up))
  
  case4.dn.genes <- Reduce(intersect, list(df.CvsA.sig[df.CvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] < -change_cutoff, unique_id],
                                           df.BvsA.nsig[, unique_id]))
  
  case4.dn <- as.data.frame(df.CvsA.sig[case4.dn.genes, ])
  case4.dn$trend <- rep("dn", nrow(case4.dn))
  
  ###########################################
  # Case 5: Significant in data 1 and data 3
  ###########################################
  # These genes would be consider markers of the first condition since they are
  # dysregulated compared to the baseline and to the second condition only
  
  case5.up.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] < -change_cutoff, unique_id],
                                           df.CvsA.nsig[, unique_id]))
  
  case5.up <- as.data.frame(df.BvsA.sig[case5.up.genes, ])
  case5.up$trend <- rep("up", nrow(case5.up))
  
  case5.dn.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] > change_cutoff, unique_id],
                                           df.CvsA.nsig[, unique_id]))
  
  case5.dn <- as.data.frame(df.BvsA.sig[case5.dn.genes, ])
  case5.dn$trend <- rep("dn", nrow(case5.dn))
  
  ###########################################
  # Case 6: Significant in data 1 and data 2
  ###########################################
  # According to the study's design, lets consider that genes dysregulated
  # in both conditions when compared to the baseline are the ones to focus on
  
  case6.up.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] > change_cutoff, unique_id],
                                           df.BvsC.nsig[, unique_id]))
  
  case6.up <- as.data.frame(df.BvsA.sig[case6.up.genes, ])
  case6.up$trend <- rep("up", nrow(case6.up))
  
  case6.dn.genes <- Reduce(intersect, list(df.BvsA.sig[df.BvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] < -change_cutoff, unique_id],
                                           df.BvsC.nsig[, unique_id]))
  
  case6.dn <- as.data.frame(df.BvsA.sig[case6.dn.genes, ])
  case6.dn$trend <- rep("dn", nrow(case6.dn))
  
  
  ################ Significant in only one comparison or neither ################
  
  # None of these cases are insightful or provide relevant information for the
  # analysis performed
  
  ########################################
  # Case 7: Significant in data 3 only
  ########################################
  
  case7.up.genes <- Reduce(intersect, list(df.BvsA.nsig[, unique_id],
                                           df.CvsA.nsig[, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] > change_cutoff, unique_id]))
  
  case7.up <- as.data.frame(df.BvsC.sig[case7.up.genes, ])
  case7.up$trend <- rep("up", nrow(case7.up))
  
  case7.dn.genes <- Reduce(intersect, list(df.BvsA.nsig[, unique_id],
                                           df.CvsA.nsig[, unique_id],
                                           df.BvsC.sig[df.BvsC.sig[, change_var] < -change_cutoff, unique_id]))
  
  case7.dn <- as.data.frame(df.BvsC.sig[case7.dn.genes, ])
  case7.dn$trend <- rep("dn", nrow(case7.dn))
  
  ########################################
  # Case 8: Significant in data 2 only
  ########################################
  
  case8.up.genes <- Reduce(intersect, list(df.BvsA.nsig[, unique_id],
                                           df.BvsC.nsig[, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] > change_cutoff, unique_id]))
  
  case8.up <- as.data.frame(df.CvsA.sig[case8.up.genes, ])
  case8.up$trend <- rep("up", nrow(case8.up))
  
  case8.dn.genes <- Reduce(intersect, list(df.BvsA.nsig[, unique_id],
                                           df.BvsC.nsig[, unique_id],
                                           df.CvsA.sig[df.CvsA.sig[, change_var] < -change_cutoff, unique_id]))
  
  case8.dn <- as.data.frame(df.CvsA.sig[case8.dn.genes, ])
  case8.dn$trend <- rep("dn", nrow(case8.dn))
  
  ########################################
  # Case 9: Significant in data 1 only
  ########################################
  
  case9.up.genes <- Reduce(intersect, list(df.CvsA.nsig[, unique_id],
                                           df.BvsC.nsig[, unique_id],
                                           df.BvsA.sig[df.BvsA.sig[, change_var] > change_cutoff, unique_id]))
  
  case9.up <- as.data.frame(df.BvsA.sig[case9.up.genes, ])
  case9.up$trend <- rep("up", nrow(case9.up))
  
  case9.dn.genes <- Reduce(intersect, list(df.CvsA.nsig[, unique_id],
                                           df.BvsC.nsig[, unique_id],
                                           df.BvsA.sig[df.BvsA.sig[, change_var] < -change_cutoff, unique_id]))
  
  case9.dn <- as.data.frame(df.BvsA.sig[case9.dn.genes, ])
  case9.dn$trend <- rep("dn", nrow(case9.dn))
  
  ########################################
  # Case 10: Significant in none
  ########################################
  
  case10.genes <- Reduce(intersect, list(df.BvsA.nsig[, unique_id],
                                         df.CvsA.nsig[, unique_id],
                                         df.BvsC.nsig[, unique_id]))
  
  case10 <- as.data.frame(df.BvsA.nsig[case10.genes, ])
  
  # Create data frames of results per cases
  case1 <- rbind(case1.up, case1.dn)
  case2 <- rbind(case2.up, case2.dn)
  case3 <- rbind(case3.up, case3.dn)
  case4 <- rbind(case4.up, case4.dn)
  case5 <- rbind(case5.up, case5.dn)
  case6 <- rbind(case6.up, case6.dn)
  case7 <- rbind(case7.up, case7.dn)
  case8 <- rbind(case8.up, case8.dn)
  case9 <- rbind(case9.up, case9.dn)
  
  return(list(Case1 = case1, Case2 = case2, Case3 = case3, Case4 = case4, Case5 = case5,
              Case6 = case6, Case7 = case7, Case8 = case8, Case9 = case9, Case10 = case10))
  
}
