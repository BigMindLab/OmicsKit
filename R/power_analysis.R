###########################
# Function power_analysis #
###########################

#' Power analysis for RNA-seq differential expression with optional plotting
#'
#' @param effect_size Numeric. Log2 fold-change to detect.
#' @param dispersion Numeric. Biological variance (σ²).
#' @param n_genes Integer. Total number of genes tested.
#' @param prop_de Numeric. Proportion of truly DE genes (0–1).
#' @param alpha Numeric. Desired FDR (type I error rate).
#' @param power_target Numeric. Desired statistical power (1 – β).
#' @param max_n Integer. Maximum sample size per group to explore.
#' @param plot Logical. If TRUE, draws the power curve; if FALSE, skips plotting.
#' @import ggplot2
#' @importFrom rlang .data
#' @export

power_analysis <- function(
    effect_size = 1,                 # log2 fold change
    dispersion = 0.1,                # biological coefficient of variation squared
    n_genes = 20000,                 # total number of genes
    prop_de = 0.05,                  # proportion of DE genes
    alpha = 0.05,                    # desired FDR
    power_target = 0.8,              # desired power
    max_n = 20,                      # maximum sample size per group
    plot = TRUE)

{
  # Calculate variance from dispersion
  sigma_sq <- dispersion

  # Z-values for significance level (Bonferroni correction as approximation for FDR control)
  m1 <- n_genes * prop_de
  m0 <- n_genes - m1
  z_alpha <- stats::qnorm(1 - alpha / 2 / m0) # adjusted for multiple testing
  z_beta <- stats::qnorm(power_target)

  # Function to compute power analytically per sample size
  compute_power <- function(n) {

    # Total variance: sigma^2 / n per group (2 groups: control vs treatment)
    se <- sqrt(2 * sigma_sq / n)
    z_effect <- abs(effect_size) / se
    power <- stats::pnorm(z_effect - z_alpha) + stats::pnorm(-z_effect - z_alpha)
    return(power)
  }

  # Compute power for range of sample sizes
  sample_sizes <- 2:max_n
  powers <- sapply(sample_sizes, compute_power)
  df <- data.frame(SampleSize = sample_sizes, Power = powers)

  # Find minimum sample size needed to reach desired power
  min_n_required <- min(df$SampleSize[df$Power >= power_target], na.rm = TRUE)

  # Plot
  if (plot) {

    p <- ggplot(df, aes(x = .data[["SampleSize"]], y = .data[["Power"]])) +
      geom_line(linewidth = 1.2, color = "#2c3e50") +
      geom_hline(yintercept = power_target, linetype = "dashed", color = "red") +
      geom_vline(xintercept = min_n_required, linetype = "dashed", color = "blue") +
      labs(title = "Power Analysis for RNA-seq Differential Expression",
           subtitle = paste0("Effect Size = ", effect_size,
                             ", Dispersion = ", dispersion,
                             ", FDR = ", alpha,
                             ", DE genes: ", prop_de * 100, "%"),
           x = "Sample Size per Group", y = "Statistical Power") +
      theme_bw(base_size = 14)

    print(p)
  }

  return(list(min_sample_size = min_n_required, power_table = df))
}
