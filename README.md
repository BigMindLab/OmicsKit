
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OmicsKit

<!-- badges: start -->
<!-- badges: end -->

The goal of OmicsKit is to help in manipulating tables and generating
plots for multi-omics analyses including genomics, transcriptomics,
proteomics, methylomics and immunoinformatics.

## Installation

You can install the development version of OmicsKit from
[GitHub](https://github.com/) with:

``` r
# Install the remotes package if needed
install.packages("remotes")

# Install from GitHub
remotes::install_github("BigMindLab/OmicsKit")

# Call library for usage
library(OmicsKit)
```

## Key features

- **Gene Annotation**: Retrieve information from Ensembl and BioMart to
  annotate gene counts tables, including transcript and gene names,
  genomic coordinates, and cross-references from various annotation
  databases.

``` r
# Example on generating transcript annotations file from Ensembl release 103

tx2gene <- get_annotations(rownames(txi$counts),
                           version = "103",
                           format = "xlsx")
```

    #>          transcriptID          geneID  symbol        biotype chromosome
    #> 1  ENST00000000233.10 ENSG00000004059    ARF5 protein_coding          7
    #> 2   ENST00000000412.8 ENSG00000003056    M6PR protein_coding         12
    #> 3  ENST00000000442.11 ENSG00000173153   ESRRA protein_coding         11
    #> 4   ENST00000001008.6 ENSG00000004478   FKBP4 protein_coding         12
    #> 5   ENST00000001146.7 ENSG00000003137 CYP26B1 protein_coding          2
    #> 6   ENST00000002125.9 ENSG00000003509 NDUFAF7 protein_coding          2
    #> 7  ENST00000002165.11 ENSG00000001036   FUCA2 protein_coding          6
    #> 8  ENST00000002501.11 ENSG00000003249  DBNDD1 protein_coding         16
    #> 9   ENST00000002596.6 ENSG00000002587  HS3ST1 protein_coding          4
    #> 10  ENST00000002829.8 ENSG00000001617  SEMA3F protein_coding          3
    #>    gene_start  gene_end gene_length
    #> 1   127588386 127591700        3315
    #> 2     8940361   8949761        9401
    #> 3    64305497  64316743       11247
    #> 4     2794970   2805423       10454
    #> 5    72129238  72147862       18625
    #> 6    37231631  37253403       21773
    #> 7   143494812 143511720       16909
    #> 8    90004871  90019890       15020
    #> 9    11393150  11429564       36415
    #> 10   50155045  50189075       34031
    #>                                                                                     description
    #> 1                                   ADP ribosylation factor 5 [Source:HGNC Symbol;Acc:HGNC:658]
    #> 2             mannose-6-phosphate receptor, cation dependent [Source:HGNC Symbol;Acc:HGNC:6752]
    #> 3                            estrogen related receptor alpha [Source:HGNC Symbol;Acc:HGNC:3471]
    #> 4                                    FKBP prolyl isomerase 4 [Source:HGNC Symbol;Acc:HGNC:3720]
    #> 5            cytochrome P450 family 26 subfamily B member 1 [Source:HGNC Symbol;Acc:HGNC:20581]
    #> 6  NADH:ubiquinone oxidoreductase complex assembly factor 7 [Source:HGNC Symbol;Acc:HGNC:28816]
    #> 7                                       alpha-L-fucosidase 2 [Source:HGNC Symbol;Acc:HGNC:4008]
    #> 8                             dysbindin domain containing 1 [Source:HGNC Symbol;Acc:HGNC:28455]
    #> 9           heparan sulfate-glucosamine 3-sulfotransferase 1 [Source:HGNC Symbol;Acc:HGNC:5194]
    #> 10                                            semaphorin 3F [Source:HGNC Symbol;Acc:HGNC:10728]

- **Dimensionality Reduction**: Generate a range of visually appealing
  plots for high-dimensional data. Includes unsupervised clustering
  methods as well.
  - PCA (Principal Component Analysis)

``` r
nice_PCA(object = transf.data,
         PCs = c(1, 2),
         ntop =  200,
         returnData = FALSE,
         variables = c(fill = "group", shape = "sex"),
         legend_names = c(fill = "Group", shape = "Sex"),
         size = 9,
         alpha = 1,
         shapes = 21:22,
         colors = my_colors,
         legend_title = 10,
         legend_elements = 8,
         legend_pos = c(0.80, 0.80, "right"),
         labels = c(var = "id", size = 3))
```

<img src="man/figures/README-pca-1.png" width="80%" /> + tSNE
(t-distributed Stochastic Neighbor Embedding)

``` r
nice_tSNE(object = transf.data,
          seed = 0,
          perplexity = 3,
          max_iterations = 10000,
          returnData = FALSE,
          variables = c(fill = "group", shape = "sex"),
          legend_names = c(fill = "Group", shape = "Sex"),
          size = 9,
          alpha = 1,
          shapes = 21:22,
          colors = my_colors,
          legend_title = 10,
          legend_elements = 8,
          legend_pos = c(0.80, 0.80, "right"),
          labels = c(var = "num", size = 3))
```

<img src="man/figures/README-tsne-1.png" width="80%" /> + UMAP (Uniform
Manifold Approximation and Projection)

``` r
nice_UMAP(object = transf.data,
          neighbors = 4,
          components = 3,
          epochs = 10000,
          returnData = FALSE,
          variables = c(fill = "group", shape = "sex"),
          legend_names = c(fill = "Group", shape = "Sex"),
          shapes = 21:22,
          colors = my_colors,
          size = 9,
          alpha = 1,
          legend_title = 10,
          legend_elements = 8,
          legend_pos = c(0.80, 0.80, "right"),
          labels = c(var = "num", size = 3))
```

<img src="man/figures/README-umap-1.png" width="80%" />

- **Counts Normalization**: Compute and extract normalized counts such
  as TPM, RPKM, FPKM, and the normalized counts from DESeq2.

``` r
# Retrieve TPMs
gene.tpm <- tpm(raw_counts = counts.gene,
                gene_lengths = counts.gene_annotations$gene_length)

# Convert to data frame
gene.tpm <- data.frame(gene.tpm)

# Add annotations
gene.tpm.annotated <- add_annotations(object = gene.tpm,
                                      reference = geneID.details,
                                      variables = annotations)
```

    #>                        S7505      S7588       S7644          geneID   symbol
    #> ENSG00000000003  0.007897175 0.01332420 0.005832601 ENSG00000000003   TSPAN6
    #> ENSG00000000005  0.000000000 0.00000000 0.000000000 ENSG00000000005     TNMD
    #> ENSG00000000419  0.425217339 0.31342321 0.266468292 ENSG00000000419     DPM1
    #> ENSG00000000457  0.319121112 0.22306157 0.152358455 ENSG00000000457    SCYL3
    #> ENSG00000000460  0.008210803 0.00826710 0.003325549 ENSG00000000460 C1orf112
    #> ENSG00000000938 11.005523543 5.41430610 4.858797243 ENSG00000000938      FGR
    #> ENSG00000000971  0.018119278 0.03191976 0.023222259 ENSG00000000971      CFH
    #> ENSG00000001036  0.333962388 0.22842915 0.191101243 ENSG00000001036    FUCA2
    #> ENSG00000001084  0.096533971 0.11264744 0.085116190 ENSG00000001084     GCLC
    #> ENSG00000001167  0.617288475 0.58584127 0.330812387 ENSG00000001167     NFYA
    #>                        biotype chromosome gene_start  gene_end gene_length
    #> ENSG00000000003 protein_coding          X  100627108 100639991       12884
    #> ENSG00000000005 protein_coding          X  100584936 100599885       14950
    #> ENSG00000000419 protein_coding         20   50934867  50958555       23689
    #> ENSG00000000457 protein_coding          1  169849631 169894267       44637
    #> ENSG00000000460 protein_coding          1  169662007 169854080      192074
    #> ENSG00000000938 protein_coding          1   27612064  27635185       23122
    #> ENSG00000000971 protein_coding          1  196652043 196747504       95462
    #> ENSG00000001036 protein_coding          6  143494812 143511720       16909
    #> ENSG00000001084 protein_coding          6   53497341  53616970      119630
    #> ENSG00000001167 protein_coding          6   41072945  41099976       27032
    #>                                                                                                    description
    #> ENSG00000000003                                              tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
    #> ENSG00000000005                                                tenomodulin [Source:HGNC Symbol;Acc:HGNC:17757]
    #> ENSG00000000419 dolichyl-phosphate mannosyltransferase subunit 1, catalytic [Source:HGNC Symbol;Acc:HGNC:3005]
    #> ENSG00000000457                                   SCY1 like pseudokinase 3 [Source:HGNC Symbol;Acc:HGNC:19285]
    #> ENSG00000000460                        chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
    #> ENSG00000000938              FGR proto-oncogene, Src family tyrosine kinase [Source:HGNC Symbol;Acc:HGNC:3697]
    #> ENSG00000000971                                         complement factor H [Source:HGNC Symbol;Acc:HGNC:4883]
    #> ENSG00000001036                                        alpha-L-fucosidase 2 [Source:HGNC Symbol;Acc:HGNC:4008]
    #> ENSG00000001084                 glutamate-cysteine ligase catalytic subunit [Source:HGNC Symbol;Acc:HGNC:4311]
    #> ENSG00000001167                nuclear transcription factor Y subunit alpha [Source:HGNC Symbol;Acc:HGNC:7804]

- **Differential Expression Results**: Filter and export differential
  expression analysis results into MS Excel or CSV formats. Filtering
  criteria include:
  - Expression change (log2 fold change).
  - Significance (False Discovery Rate, FDR).
  - Detectability (`Requena et al., 2024, Nat. Comms.`).

``` r
detect_list <- detect_filter(norm.counts = normalized.counts[, 1:21],
                             df.BvsA = res.T_N.sig,
                             df.CvsA = res.M_N.sig,
                             df.DvsA = NULL,
                             cutoffs = c(50, 50, 0),
                             samples.baseline = 1:3,
                             samples.condition1 = 4:6,
                             samples.condition2 = 7:9,
                             samples.condition3 = NULL)
```

    #> $Comparison1
    #>                  baseMean log2FoldChange     lfcSE      stat       pvalue
    #> ENSG00000015532 257.21417     -1.5441220 0.4673795 -3.303786 9.538854e-04
    #> ENSG00000040633 439.99903      0.4652085 0.1393893  3.337477 8.454280e-04
    #> ENSG00000043514  74.19499     -0.7064732 0.2003565 -3.526082 4.217567e-04
    #> ENSG00000050748 259.94720     -0.6654422 0.1947707 -3.416541 6.342210e-04
    #> ENSG00000067840 499.56771     -1.5485951 0.3585761 -4.318735 1.569261e-05
    #> ENSG00000073417 126.05703     -1.2149724 0.2967272 -4.094577 4.229395e-05
    #>                       padj         ensembl symbol        biotype chromosome
    #> ENSG00000015532 0.24746629 ENSG00000015532  XYLT2 protein_coding         17
    #> ENSG00000040633 0.24746629 ENSG00000040633  PHF23 protein_coding         17
    #> ENSG00000043514 0.19637741 ENSG00000043514  TRIT1 protein_coding          1
    #> ENSG00000050748 0.22608000 ENSG00000050748  MAPK9 protein_coding          5
    #> ENSG00000067840 0.04826785 ENSG00000067840  PDZD4 protein_coding          X
    #> ENSG00000073417 0.08672609 ENSG00000073417  PDE8A protein_coding         15
    #>                 gene_start  gene_end gene_length
    #> ENSG00000015532   50346126  50363138       17013
    #> ENSG00000040633    7235029   7239722        4694
    #> ENSG00000043514   39838110  39883511       45402
    #> ENSG00000050748  180233143 180292099       58957
    #> ENSG00000067840  153802166 153830565       28400
    #> ENSG00000073417   84980440  85139145      158706
    #>                                                                           description
    #> ENSG00000015532              xylosyltransferase 2 [Source:HGNC Symbol;Acc:HGNC:15517]
    #> ENSG00000040633             PHD finger protein 23 [Source:HGNC Symbol;Acc:HGNC:28428]
    #> ENSG00000043514     tRNA isopentenyltransferase 1 [Source:HGNC Symbol;Acc:HGNC:20286]
    #> ENSG00000050748 mitogen-activated protein kinase 9 [Source:HGNC Symbol;Acc:HGNC:6886]
    #> ENSG00000067840           PDZ domain containing 4 [Source:HGNC Symbol;Acc:HGNC:21167]
    #> ENSG00000073417               phosphodiesterase 8A [Source:HGNC Symbol;Acc:HGNC:8793]
    #> 
    #> $DetectGenes
    #>  [1] "ENSG00000015532" "ENSG00000040633" "ENSG00000043514" "ENSG00000050748"
    #>  [5] "ENSG00000067840" "ENSG00000073417" "ENSG00000092820" "ENSG00000102317"
    #>  [9] "ENSG00000109519" "ENSG00000115289" "ENSG00000115607" "ENSG00000118503"
    #> [13] "ENSG00000118507" "ENSG00000118520" "ENSG00000121057" "ENSG00000123329"
    #> [17] "ENSG00000125835" "ENSG00000134256" "ENSG00000138777" "ENSG00000145241"
    #> [21] "ENSG00000145675" "ENSG00000146859" "ENSG00000147459" "ENSG00000148488"
    #> [25] "ENSG00000152061" "ENSG00000154640" "ENSG00000155749" "ENSG00000156675"
    #> [29] "ENSG00000160633" "ENSG00000164284" "ENSG00000167118" "ENSG00000168887"
    #> [33] "ENSG00000169136" "ENSG00000171865" "ENSG00000174007" "ENSG00000181026"
    #> [37] "ENSG00000181523" "ENSG00000183484" "ENSG00000183726" "ENSG00000186407"
    #> [41] "ENSG00000188315" "ENSG00000197448" "ENSG00000198841" "ENSG00000234389"
    #> [45] "ENSG00000267283" "ENSG00000078589" "ENSG00000100416" "ENSG00000101191"
    #> [49] "ENSG00000135164" "ENSG00000135763"
    #> 
    #> $Comparison2
    #>                  baseMean log2FoldChange     lfcSE      stat       pvalue
    #> ENSG00000078589 456.36884     -0.5806294 0.1631107 -3.559726 3.712419e-04
    #> ENSG00000100416 191.21091     -0.8072890 0.1874153 -4.307487 1.651198e-05
    #> ENSG00000101191 807.06214     -0.4319645 0.1163969 -3.711135 2.063322e-04
    #> ENSG00000121057 161.09141     -1.0596203 0.2722106 -3.892649 9.915571e-05
    #> ENSG00000135164 713.36398     -2.4825790 0.6778906 -3.662212 2.500468e-04
    #> ENSG00000135763  55.37574     -1.1755918 0.3069112 -3.830397 1.279366e-04
    #>                       padj         ensembl symbol        biotype chromosome
    #> ENSG00000078589 0.24733759 ENSG00000078589 P2RY10 protein_coding          X
    #> ENSG00000100416 0.03650248 ENSG00000100416   TRMU protein_coding         22
    #> ENSG00000101191 0.19621498 ENSG00000101191  DIDO1 protein_coding         20
    #> ENSG00000121057 0.15908922 ENSG00000121057  AKAP1 protein_coding         17
    #> ENSG00000135164 0.21630136 ENSG00000135164  DMTF1 protein_coding          7
    #> ENSG00000135763 0.15908922 ENSG00000135763   URB2 protein_coding          1
    #>                 gene_start  gene_end gene_length
    #> ENSG00000078589   78945332  78963727       18396
    #> ENSG00000100416   46330875  46357340       26466
    #> ENSG00000101191   62877738  62937952       60215
    #> ENSG00000121057   57085092  57121346       36255
    #> ENSG00000135164   87152409  87196337       43929
    #> ENSG00000135763  229626247 229660200       33954
    #>                                                                                          description
    #> ENSG00000078589                    P2Y receptor family member 10 [Source:HGNC Symbol;Acc:HGNC:19906]
    #> ENSG00000100416               tRNA mitochondrial 2-thiouridylase [Source:HGNC Symbol;Acc:HGNC:25481]
    #> ENSG00000101191                       death inducer-obliterator 1 [Source:HGNC Symbol;Acc:HGNC:2680]
    #> ENSG00000121057                       A-kinase anchoring protein 1 [Source:HGNC Symbol;Acc:HGNC:367]
    #> ENSG00000135164 cyclin D binding myb like transcription factor 1 [Source:HGNC Symbol;Acc:HGNC:14603]
    #> ENSG00000135763                 URB2 ribosome biogenesis homolog [Source:HGNC Symbol;Acc:HGNC:28967]
    #> 
    #> $Comparison3
    #>                  baseMean log2FoldChange     lfcSE      stat       pvalue
    #> ENSG00000004777  59.42425     -0.7664742 0.1664979 -4.603506 4.154375e-06
    #> ENSG00000014914  50.00342      0.9984404 0.2941904  3.393858 6.891546e-04
    #> ENSG00000015532 257.21417     -1.5853106 0.3652982 -4.339771 1.426310e-05
    #> ENSG00000040633 439.99903      0.4286551 0.1094576  3.916175 8.996484e-05
    #> ENSG00000051108 611.57428     -0.5931505 0.1560919 -3.800009 1.446909e-04
    #> ENSG00000070366 309.66550     -0.5852590 0.1641947 -3.564421 3.646608e-04
    #>                        padj         ensembl   symbol        biotype chromosome
    #> ENSG00000004777 0.009346896 ENSG00000004777 ARHGAP33 protein_coding         19
    #> ENSG00000014914 0.229188351 ENSG00000014914   MTMR11 protein_coding          1
    #> ENSG00000015532 0.019606051 ENSG00000015532    XYLT2 protein_coding         17
    #> ENSG00000040633 0.077291042 ENSG00000040633    PHF23 protein_coding         17
    #> ENSG00000051108 0.114745426 ENSG00000051108  HERPUD1 protein_coding         16
    #> ENSG00000070366 0.167326005 ENSG00000070366     SMG6 protein_coding         17
    #>                 gene_start  gene_end gene_length
    #> ENSG00000004777   35774532  35788822       14291
    #> ENSG00000014914  149928651 149936879        8229
    #> ENSG00000015532   50346126  50363138       17013
    #> ENSG00000040633    7235029   7239722        4694
    #> ENSG00000051108   56932142  56944864       12723
    #> ENSG00000070366    2059839   2303785      243947
    #>                                                                                                        description
    #> ENSG00000004777                               Rho GTPase activating protein 33 [Source:HGNC Symbol;Acc:HGNC:23085]
    #> ENSG00000014914                                myotubularin related protein 11 [Source:HGNC Symbol;Acc:HGNC:24307]
    #> ENSG00000015532                                           xylosyltransferase 2 [Source:HGNC Symbol;Acc:HGNC:15517]
    #> ENSG00000040633                                          PHD finger protein 23 [Source:HGNC Symbol;Acc:HGNC:28428]
    #> ENSG00000051108 homocysteine inducible ER protein with ubiquitin like domain 1 [Source:HGNC Symbol;Acc:HGNC:13744]
    #> ENSG00000070366                       SMG6 nonsense mediated mRNA decay factor [Source:HGNC Symbol;Acc:HGNC:17809]

- **Case Organization**: Automatically categorize results from three
  pairwise differential expression analyses or Gene Set Enrichment
  Analysis (e.g., B vs A, C vs A, C vs B) into 10 mutually exclusive
  cases (`BigMind, 2024, manuscript in preparation`).

``` r
DEGs_sig <- split_cases(df.BvsA = res.T_N,
                        df.CvsA = res.M_N,
                        df.BvsC = res.M_T,
                        unique_id = "ensembl",
                        significance_var = "padj",
                        significance_cutoff = 0.25,
                        change_var = "log2FoldChange",
                        change_cutoff = 0)

# Filter the whole detectability list by a new threshold
for (i in names(DEGs_sig)) {
  DEGs_sig[[i]] <- DEGs_sig[[i]][rownames(DEGs_sig[[i]]) %in% detect_list$DetectGenes, ]
  DEGs_sig[[i]] <- DEGs_sig[[i]][DEGs_sig[[i]]$padj < 0.05, ]
}
```

    #> $Case6
    #>                 baseMean log2FoldChange     lfcSE      stat       pvalue
    #> ENSG00000109519 158.4705     -0.7058074 0.1620767 -4.354774 1.332044e-05
    #>                       padj         ensembl symbol        biotype chromosome
    #> ENSG00000109519 0.04826785 ENSG00000109519 GRPEL1 protein_coding          4
    #>                 gene_start gene_end gene_length
    #> ENSG00000109519    7058895  7068064        9170
    #>                                                                    description
    #> ENSG00000109519 GrpE like 1, mitochondrial [Source:HGNC Symbol;Acc:HGNC:19696]
    #>                 trend
    #> ENSG00000109519    dn
    #> 
    #> $Case8
    #>                  baseMean log2FoldChange     lfcSE      stat       pvalue
    #> ENSG00000130656 124.20885      4.3999189 0.9629323  4.569292 4.893748e-06
    #> ENSG00000004777  59.42425     -0.7664742 0.1664979 -4.603506 4.154375e-06
    #> ENSG00000137842  97.42166     -0.6824943 0.1504292 -4.536981 5.706537e-06
    #>                        padj         ensembl   symbol        biotype chromosome
    #> ENSG00000130656 0.009346896 ENSG00000130656      HBZ protein_coding         16
    #> ENSG00000004777 0.009346896 ENSG00000004777 ARHGAP33 protein_coding         19
    #> ENSG00000137842 0.009346896 ENSG00000137842   TMEM62 protein_coding         15
    #>                 gene_start gene_end gene_length
    #> ENSG00000130656     142728   154503       11776
    #> ENSG00000004777   35774532 35788822       14291
    #> ENSG00000137842   43123279 43185144       61866
    #>                                                                          description
    #> ENSG00000130656           hemoglobin subunit zeta [Source:HGNC Symbol;Acc:HGNC:4835]
    #> ENSG00000004777 Rho GTPase activating protein 33 [Source:HGNC Symbol;Acc:HGNC:23085]
    #> ENSG00000137842         transmembrane protein 62 [Source:HGNC Symbol;Acc:HGNC:26269]
    #>                 trend
    #> ENSG00000130656    up
    #> ENSG00000004777    dn
    #> ENSG00000137842    dn
    #> 
    #> $Case9
    #>                  baseMean log2FoldChange     lfcSE      stat       pvalue
    #> ENSG00000067840  499.5677     -1.5485951 0.3585761 -4.318735 1.569261e-05
    #> ENSG00000102317 1755.7926     -0.5060843 0.1169605 -4.326966 1.511769e-05
    #>                       padj         ensembl symbol        biotype chromosome
    #> ENSG00000067840 0.04826785 ENSG00000067840  PDZD4 protein_coding          X
    #> ENSG00000102317 0.04826785 ENSG00000102317   RBM3 protein_coding          X
    #>                 gene_start  gene_end gene_length
    #> ENSG00000067840  153802166 153830565       28400
    #> ENSG00000102317   48574449  48581162        6714
    #>                                                                    description
    #> ENSG00000067840    PDZ domain containing 4 [Source:HGNC Symbol;Acc:HGNC:21167]
    #> ENSG00000102317 RNA binding motif protein 3 [Source:HGNC Symbol;Acc:HGNC:9900]
    #>                 trend
    #> ENSG00000067840    dn
    #> ENSG00000102317    dn

- **Common plots**: Generate a range of visually appealing plots to
  display differential expression analysis results. Here are some
  examples
  - **Volcano plots**

<img src="man/figures/README-Volcano_plot.jpg" width="100%" /> +
**Heatmaps**

``` r
heatmap_Ap_An
#> $tree_row
#> 
#> Call:
#> hclust(d = d, method = method)
#> 
#> Cluster method   : mcquitty 
#> Distance         : euclidean 
#> Number of objects: 24 
#> 
#> 
#> $tree_col
#> [1] NA
#> 
#> $kmeans
#> [1] NA
#> 
#> $gtable
#> TableGrob (5 x 6) "layout": 8 grobs
#>   z     cells                 name                          grob
#> 1 1 (4-4,1-1)             row_tree polyline[GRID.polyline.22146]
#> 2 2 (4-4,3-3)               matrix       gTree[GRID.gTree.22148]
#> 3 3 (5-5,3-3)            col_names         text[GRID.text.22149]
#> 4 4 (4-4,4-4)            row_names         text[GRID.text.22150]
#> 5 5 (3-3,3-3)       col_annotation         rect[GRID.rect.22151]
#> 6 6 (3-3,4-4) col_annotation_names         text[GRID.text.22152]
#> 7 7 (3-5,6-6)    annotation_legend       gTree[GRID.gTree.22157]
#> 8 8 (3-5,5-5)               legend       gTree[GRID.gTree.22160]
#> 
#> attr(,"class")
#> [1] "pheatmap"
```

    + **Enrichment plots**

<img src="man/figures/README-Balloon_plot.jpeg" width="100%" />

    + **Box-Scatter-Violin (BSV) plots**

<img src="man/figures/README-BSV_plot.jpeg" width="100%" />

## Examples

For a more detailed workflow on Differential Expression Analysis with
the application of the `OmicsKit` suit please check the custom
[BigMind](https://github.com/BigMindLab)’s pipeline for
[DESeq2](https://github.com/BigMindLab/DESeq2).

## Developed by

- David R. Requena Anicama, Ph.D.

  - Author’s name: David Requena
  - [Google
    Scholar](https://scholar.google.com/citations?user=uI01iS4AAAAJ&hl=en)
  - [ORCID: 0000-0002-5968-1133](https://orcid.org/0000-0002-5968-1133)

- Daniel F. Guevara Díaz, B.Sc.(s)

  - Author’s name: Daniel F. Guevara-Díaz
  - [Google
    Scholar](https://scholar.google.com/citations?hl=en&user=tqT7vr8AAAAJ)
  - [ORCID: 0009-0001-2786-8729](https://orcid.org/0009-0001-2786-8729)

## License

CC BY-NC-SA 4.0

## Contact

<david.requena@nyulangone.org>
