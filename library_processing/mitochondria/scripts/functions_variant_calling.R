# Here I modify the generation of a ternary matrix. In the original mitoClone package
# only the allele frequency was considered when assigning cells as reference and 
# mutant. Here I add a minimum threshold coverage, so that cells with less than a certain 
# number of reads are labelled as dropouts


# add coverage threshold to mutationCalls function. The default will be 10 reads
mutationCallsFromMatrixCoverage <- function (M, N, cluster = NULL, metadata = data.frame(row.names = rownames(M)), 
                                             binarize = 0.05, min_coverage = 10) 
{
  colnames(M) <- make.names(colnames(M))
  colnames(N) <- make.names(colnames(N))
  binfun <- function(M, N) {
    alleleRatio <- M/(M + N)
    coverage <- M+N
    ternary_list <- lapply(1:ncol(M), function(x){
      
      data <- data.frame(af = alleleRatio[,x], cov_site = coverage[,x]) %>% 
                mutate(ternary = ifelse(cov_site < min_coverage, "?", ifelse(af>binarize, "1", "0"))) %>% 
                dplyr::select(ternary)
    })
    
    ternary_matrix <- do.call("cbind", ternary_list)
    colnames(ternary_matrix) <- colnames(M)
    
    return(as.matrix(ternary_matrix))
    
}
  out <- new("mutationCalls", M = M, N = N, metadata = metadata, 
             ternary = binfun(M, N))
  if (!is.null(cluster)) 
    out@cluster <- cluster
  else {
    out@cluster <- apply(out@ternary != "?", 2, mean) > 
      0.2
  }
  out
}


mutationCallsFromBlacklistCoverage <- function (BaseCounts, lim.cov = 20, min.af = 0.2, min.num.samples = 0.01 * 
                                                  length(BaseCounts), min.af.universal = min.af, universal.var.cells = 0.95 * 
                                                  length(BaseCounts), blacklists.use = blacklists, max.var.na = 0.5, 
                                                max.cell.na = 0.95, ...) 
{
  varaf <- parallel::mclapply(BaseCounts, function(x) {
    x <- x[, 1:4]
    zeroes <- rowSums(x) < lim.cov
    x.af <- x/(x + apply(x, 1, max))
    x.af <- reshape2::melt(x.af)
    colnames(x.af) <- c("pos", "nt", "af")
    x.af <- x.af[!(mito.dna[x.af$pos] == x.af$nt), ]
    x.af <- x.af[!(mito.dna[x.af$pos] == "N"), ]
    x.af$name <- paste0(x.af$pos, " ", mito.dna[x.af$pos], 
                        ">", x.af$nt)
    x.af$af[x.af$pos %in% which(zeroes)] <- NA
    x <- x.af$af
    names(x) <- x.af$name
    return(x)
  }, mc.cores = 10)
  varaf <- do.call(cbind, varaf)
  varaf <- varaf[rowSums(varaf > min.af, na.rm = TRUE) >= 
                   min.num.samples, ]
  is.names <- sapply(blacklists.use, function(x) typeof(x) == 
                       "character")
  if (sum(is.names) > 0) {
    removal.names.list <- unique(unlist(blacklists.use[is.names]))
    varaf <- varaf[!row.names(varaf) %in% removal.names.list, 
                   ]
  }
  if (sum(!is.names) > 0) {
    removal.ranges.list <- unique(unlist(GenomicRanges::GRangesList(blacklists.use[!is.names])))
    varaf <- varaf[-c(S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2gr(row.names(varaf)), 
                                                                       removal.ranges.list))), ]
  }
  varaf <- varaf[rowSums(varaf, na.rm = T) > 0, ]
  varaf <- varaf[!rowSums(varaf >= min.af.universal, na.rm = TRUE) >= 
                   universal.var.cells, ]
  varaf <- varaf[rowSums(is.na(varaf)) < max.var.na * NCOL(varaf), 
                 ]
  varaf <- varaf[, colSums(is.na(varaf)) < max.cell.na * NROW(varaf)]
  MN <- pullcounts.vars(BaseCounts, rownames(varaf), colnames(varaf))
  mutationCallsFromMatrixCoverage(t(MN$M), t(MN$N), ...)
}


# colours for UMAPs

warriors_yellow <- "#FDBB30"
wine_red <- "#860038"
lakers_purple <- "#552583"
apple_gray <- "#7D7D7D"
jazz_blue <- "#002F6C"
bucks_green <- "#274E37"
dallas_blue <- "#2A6BAC"
bulls_red <- "#CE1141"
facebook_blue <- "#4267B2"


# function to compute association tests using Fisher's exact test
getasso <- function(m1,m2) {
  testm <- matrix( data = c(
    sum(m1 =="1" & m2 == "1"),
    sum(m1 =="1" & m2 == "0"),
    sum(m1 == "0" & m2 == "1"),
    sum(m1 == "0" & m2 == "0")
  ), ncol = 2  )
  fisher.test(testm)
}


# function to change the ternary function into: reference, mutated and dropout
rename_ternary <- function(x){
  
  sapply(1:length(x), function(i){
    
    name <- ""
    
    row <- x[i]
    
    if(row == "1"){name <- "mutant"}
    
    else if(row == "0"){name <- "reference"}
    
    else{name <-  "dropout"}
    
    return(name)
    
    
  })
  
}


# function to convert NA into "dropout"
rename_na <- function(x){
  
  sapply(1:length(x), function(i){
    
    name <- ""
    
    row <- x[i]
    
    if(is.na(row)){name <- "dropout"}
    
    else{name <-  row}
    
    return(name)
  })
}


# make UMAP with annotated clusters from Anne's patients (annotated by Lars and Michael)
make_umap_clusters_Anne <- function(umap_table, theme, sample){
  
  
  umap_plot <- ggplot(umap_table, aes(x = umap_1, y = umap_2)) +
                  geom_point(data = umap_table %>% filter(sample != "AKLA1b"), size = 0.5, color = "gray") +
                  geom_point(data = umap_table %>% filter(sample == "AKLA1b"), size = 1, aes(color = cell_type)) +
                  theme_classic() +
                  #scale_colour_manual(values = c(dallas_blue, bulls_red, "gray")) + 
                  theme(axis.title = element_text(size = 14),
                        axis.text = element_text(size = 12, color = "black"),
                        legend.title = element_blank(),
                        plot.title = element_text(size = 18, hjust = 0.5),
                        legend.text = element_text(size = 13),
                        legend.position = "bottom") +
                  ggtitle(paste0("Cell populations ", sample))
                  #guides(color = guide_legend(override.aes = list(size = 3)))
  
  
}


# make ternary matrix for cell population (to be used in association tests)
ternary_cell_type <- function(cell_types, table){
  
      # get ternary matrix for each cell type present in the sample
      matrix <- lapply(1:length(cell_types), function(i){
        
                ternary <- ifelse(table$cell_type == cell_types[i], "1", "0")
        
      })
      
      # merge all populations into 1 ternary matrix
      matrix <- do.call("cbind", matrix)
      
      
      # add cell type names as column names
      colnames(matrix) <- cell_types
      
      
      # add cell barcodes as rownames
      rownames(matrix) <- table$cell_barcode
      
      
      return(matrix)
      
  
}


# make ternary matrix for nuclear mutations and CNAs (to be used in association tests)
ternary_mutations <- function(mutations, table){
  
  # get ternary matrix for each cell type present in the sample
  matrix <- lapply(1:length(mutations), function(i){
    
    ternary <- ifelse(table %>% pull(mutations[i]) == "reference", "0", 
                      ifelse(table %>% pull(mutations[i]) == "dropout", "?", "1"))
    
  })
  
  # merge all populations into 1 ternary matrix
  matrix <- do.call("cbind", matrix)
  
  
  # add mutation names as column names
  colnames(matrix) <- mutations
  
  
  # add cell barcodes as rownames
  rownames(matrix) <- table$cell_barcode
  
  
  return(matrix)
  
  
}


# make mito count tables from mGATK output
make_mito_count_table <- function(sample, mito_path, output_path = here("data/")){
  
  
  # load mgatk output object
  mito_object <- readRDS(paste0(mito_path, "/", sample, "/mito_counts/count_table/final/", sample, ".rds"))
  
  
  # reshape matrices into an output compatible with mitoClone functions
  list_clean_count_tables <- mclapply(1:length(colnames(mito_object)), function(x){
    
    
    table <- data.frame(A = assays(mito_object)$A_counts_fw[,colnames(mito_object)[x]],
                        T = assays(mito_object)$T_counts_fw[,colnames(mito_object)[x]],
                        C = assays(mito_object)$C_counts_fw[,colnames(mito_object)[x]],
                        G = assays(mito_object)$G_counts_fw[,colnames(mito_object)[x]]) %>% 
      as.matrix()
    
  }, mc.cores = 8)
  
  
  # add cell barcodes as column names
  names(list_clean_count_tables) <- gsub("([A-Z]{16}).+$", "\\1", colnames(mito_object))
  
  
  # save the table
  saveRDS(list_clean_count_tables, file = paste0(output_path, "variant_count_table.rds"))
  
  
  # return list of count tables
  return(list_clean_count_tables)
  
  
}


# function to make association tests between mitochondrial mutations and nuclear mutations/cell_types
association_tests <- function(mitoClone_object, ternary_table){
  
  
  # get barcodes for which cell type is known and there is coverage of mitochondrial mutations
  valid_barcodes <- rownames(mitoClone_object@ternary)[which(rownames(mitoClone_object@ternary) %in% rownames(ternary_table))]
  
  
  # make Fisher tests 
  fisher_tests <- lapply(1:length(colnames(ternary_table)), function(i){
    
    
   # Association tests with T cells
   fisher_matrix <- apply(mitoClone_object@ternary[valid_barcodes, ], 2, getasso, 
                          m2 = ternary_table[valid_barcodes, colnames(ternary_table)[i]])
  
  
  
    data_frame <- data.frame(name = colnames(mitoClone_object@ternary),
                             pval = -log10(p.adjust(sapply(fisher_matrix, "[[", "p.value"), 
                                                    method = "BH")),
                             or = sapply(fisher_matrix, "[[", "estimate"))
    
  })
  
  
  names(fisher_tests) <- colnames(ternary_table)
  
  fisher_tests
  
}


# make qplot of associations tests
make_qplots <- function(fisher_output, cell_types = FALSE){
  
  
  # in case of associations to cell populations plots will be generated when at least one variant is significantly 
  # associated to a particularly population
  if(cell_types){
    
    sig_columns <- unlist(lapply(1:length(fisher_output), function(i){
      
                          table <- fisher_output[[i]]
                          
                          if(max(table$pval) > 1.5){population <- names(fisher_output[i])}
      
    }))
    
    
  }else{sig_columns <- names(fisher_output)}
  
  
  qplots <- lapply(1:length(sig_columns), function(i){
    
              qplot(x = name , y = pval, color = paste0(or <1, pval > 1),
                    data=na.omit(fisher_output[[sig_columns[i]]] %>% filter(pval > 1.5))) + 
                    coord_flip() + theme_classic()+ 
                    theme(panel.grid = element_blank(), 
                          axis.title.x = element_text(color="black", size = 14), 
                          axis.text = element_text(color="black", size = 12), 
                          axis.text.y = element_text(size = 12)) + 
                    xlab("") + ylab(paste0("Association w/ ",  sig_columns[i], "\n(-log10 p)"))+
                    scale_color_manual(values = c("TRUETRUE" = "blue",
                                                  "FALSETRUE" = "red", 
                                                  "FALSEFALSE" = "grey", 
                                                  "TRUEFALSE" = "grey"), guide=F)
    
  })
  
}
  
  

# make UMAPs with mutation label
mutation_umaps <- function(table, variants, sample_name){
  
  umap <- lapply(1:length(variants), function(i){
    
          site <- variants[i]
          
          # change name of mito variants (follow pattern XPOS.Ref.Alt.)
          
          if(str_split(site, "")[[1]][1] == "X" & !is.na(as.integer(str_split(site, "")[[1]][2]))){
            
            mito_variant = TRUE
            
            position <- gsub("^X(\\d+).*", "\\1", site)
            
            reference <- gsub("^.+\\d+([ATCG])\\..+$", "\\1",site)
            
            mutation <- gsub("^.*\\.([ATCG])$", "\\1",site)
            
          }else{mito_variant <- FALSE}
          
          
          plot <- ggplot(table, aes_string(x = "umap_1", y = "umap_2")) +
                    geom_point(data = table %>% filter(!sample %in% sample_name), 
                               color = "gray",
                               size = 0.25, alpha = 0.5) +
                    geom_point(data = table %>% filter(sample %in% sample_name), 
                               aes_string(color = site),
                               size = 1) +
                    theme_classic() +
                    scale_colour_manual(values = c(dallas_blue, bulls_red, "gray")) + 
                    theme(axis.title = element_text(size = 14),
                          axis.text = element_text(size = 12, color = "black"),
                          legend.title = element_blank(),
                          plot.title = element_text(size = 18, hjust = 0.5),
                          legend.text = element_text(size = 13),
                          legend.position = "bottom") +
                    ggtitle(ifelse(mito_variant, paste0("mt:", position, reference, ">", mutation), site))+
                    guides(color = guide_legend(override.aes = list(size = 3)))
    
  })
  
}

  

# make table with clonal information and confidence assignment

clone_table <- function(clustered_clones, umap_table){
  
  data.frame(cell_barcode = rownames(clustered_clones@M),
             clone = as.factor(apply(clustered_clones@mainClone, 1,
                                     which.max)),
             confidence = apply(clustered_clones@mainClone, 1, max)) %>% 
            left_join(umap_table) %>% 
           mutate(transparency = factor(case_when(
             confidence > 0.95 ~ "1", confidence > 0.8 ~ "0.8", 
             confidence > 0.7 ~ "0.7",
             confidence > 0.6 ~ "0.6", confidence > 0 ~ "0.5")))
  
}

# heatmap of clones showing allele frequency. Cells with coverage < 10 are depitcted in white as dropouts.
heatmap_clones_af <- function(clustered_clones){
  
  # compute de allele frequency of each site. Sites with coverage < 10 are set to NA.
  heatmap_data <- data.frame(ifelse(clustered_clones@M + clustered_clones@N >= 10,
                                    clustered_clones@M/(clustered_clones@M + clustered_clones@N),
                                    NA),
                             row.names = rownames(clustered_clones@M))
  
  
  # transpose the data to suit the pheatmap function
  heatmap_data <- t(heatmap_data[, getNodes(clustered_clones@tree)[-1]])
  
  
  # create an object to store annotations shown in the plot
  annos <- data.frame(row.names = rownames(clustered_clones@M), 
                      clustered_clones@metadata)
  
  
  # get the main clone assigned to each cell
  annos$mainClone <- as.factor(apply(clustered_clones@mainClone,1, which.max))
  
  
  # get the likelihood assignment for each cell
  annos$confidence <- apply(clustered_clones@mainClone, 1, max)
  
  
  # order cells into clones 
  heatmap_data <- heatmap_data[,order(annos$mainClone)]
  
  
  # plot heatmap
  heatmap_presentation <- pheatmap::pheatmap(heatmap_data, cluster_cols = F, cluster_rows = F, 
                                             show_colnames = F, 
                                             color = colorRampPalette(rev(c("#9B0000", "#FFD72E", "#FFD72E", "#00009B")))(100), 
                                             annotation_col = annos,
                                             na_col = "white")
  
}


# function to make heatmap of clones showing the mutational status
heatmap_clones_ternary <- function(clustered_clones){
  
  # compute de allele frequency of each site. Sites with coverage < 10 are set to NA.
  heatmap_data <- data.frame(apply(clustered_clones@ternary, 2, function(x) 
                              ifelse(x == "1", 1, ifelse(x == "?", 0, -1))),
                             row.names = rownames(clustered_clones@M))
  
  
  # transpose the data to suit the pheatmap function
  heatmap_data <- t(heatmap_data[, getNodes(clustered_clones@tree)[-1]])
  
  
  # create an object to store annotations shown in the plot
  annos <- data.frame(row.names = rownames(clustered_clones@M), 
                      clustered_clones@metadata)
  
  
  # get the main clone assigned to each cell
  annos$mainClone <- as.factor(apply(clustered_clones@mainClone,1, which.max))
  
  
  # get the likelihood assignment for each cell
  annos$confidence <- apply(clustered_clones@mainClone, 1, max)
  
  
  # order cells into clones 
  heatmap_data <- heatmap_data[,order(annos$mainClone)]
  
  
  # plot heatmap
  heatmap_presentation <- pheatmap::pheatmap(heatmap_data, cluster_cols = F, cluster_rows = F, 
                                             show_colnames = F, 
                                             color = colorRampPalette(rev(c("#9B0000", "#FFD72E", "#FFD72E", "#00009B")))(100), 
                                             annotation_col = annos)
  
}


# mutation calls from Cohort with min 10x coverage to assign a cell as mutant or reference in the final ternary matrix
mutationCallsFromCohortCoverage <- function (BaseCounts, patient, MINREADS = 5, MINCELL = 20, MINFRAC = 0.1, 
                                      MINCELLS.PATIENT = 10, MINRELATIVE.PATIENT = 0.01, MINRELATIVE.OTHER = 0.1) {
  nuc.count.per.position.array <- array(data = 0, dim = c(length(BaseCounts), 
                                                          nrow(BaseCounts[[1]]), ncol(BaseCounts[[1]])), dimnames = list(names(BaseCounts), 
                                                                                                                         paste0("X", 1:nrow(BaseCounts[[1]])), colnames(BaseCounts[[1]])))
  for (i in 1:length(BaseCounts)) nuc.count.per.position.array[i, 
                                                               , ] <- BaseCounts[[i]]
  sum.overall <- apply(nuc.count.per.position.array, c(2, 
                                                       3), sum)
  reference <- colnames(sum.overall)[apply(sum.overall, 1, 
                                           which.max)]
  mt.reads.per.cell <- apply(nuc.count.per.position.array, 
                             1, sum)
  variant_calls <- lapply(1:length(reference), function(pos) {
    support <- apply(nuc.count.per.position.array[, pos, 
                                                  ] > MINREADS, 2, sum)
    support <- support[!names(support) %in% c(reference[pos], 
                                              "N")]
    use <- names(support)[support > MINCELL]
    if (length(use) == 0) 
      NULL
    else {
      out <- matrix(data = NA, ncol = length(use), nrow = nrow(nuc.count.per.position.array[, 
                                                                                            pos, ]), dimnames = list(rownames(nuc.count.per.position.array[, 
                                                                                                                                                           pos, ]), paste0(pos, reference[pos], ">", use)))
      for (i in 1:length(use)) {
        pos_sum <- apply(nuc.count.per.position.array[, 
                                                      pos, ], 1, sum)
        condition_mut <- nuc.count.per.position.array[, 
                                                      pos, use[i]] > MINREADS & nuc.count.per.position.array[, 
                                                                                                             pos, use[i]] > MINFRAC * pos_sum
        condition_ref <- nuc.count.per.position.array[, 
                                                      pos, reference[pos]] > MINREADS & nuc.count.per.position.array[, 
                                                                                                                     pos, reference[pos]] > MINFRAC * pos_sum
        out[, i] <- ifelse(condition_mut, ifelse(condition_ref, 
                                                 "BOTH", "MUT"), ifelse(condition_ref, "WT", 
                                                                        "DROP"))
      }
      out
    }
  })
  variant_calls <- do.call(cbind, variant_calls)
  varcount.bypatient <- sapply(unique(patient), function(pa) {
    apply(variant_calls[patient == pa, ], 2, function(x) sum(x %in% 
                                                               c("BOTH", "MUT")))
  })
  patient.count <- as.vector(table(patient)[colnames(varcount.bypatient)])
  names(patient.count) <- colnames(varcount.bypatient)
  varcount.bypatient.fraction <- t(t(varcount.bypatient)/patient.count)
  filter <- apply(varcount.bypatient, 1, max) > MINCELLS.PATIENT & 
    apply(varcount.bypatient, 1, function(x) max(x)/patient.count[which.max(x)]) > 
    MINRELATIVE.PATIENT
  patientfilter <- filter & apply(varcount.bypatient.fraction, 
                                  1, function(x) sum(x > MINRELATIVE.OTHER * max(x))) == 
    1 & apply(varcount.bypatient, 1, function(x) sum(x >= 
                                                       MINCELLS.PATIENT)) == 1
  multipatient <- filter & apply(varcount.bypatient.fraction, 
                                 1, function(x) sum(x > MINRELATIVE.OTHER * max(x))) > 
    1 & apply(varcount.bypatient, 1, function(x) sum(x >= 
                                                       MINCELLS.PATIENT)) > 1 & !grepl(">-", rownames(varcount.bypatient))
  mutation.bypatient <- colnames(varcount.bypatient)[apply(varcount.bypatient[patientfilter, 
                                                                              ], 1, which.max)]
  variant_calls_selected <- variant_calls[, patientfilter]
  mutation.bypatient <- mutation.bypatient[!grepl("->", colnames(variant_calls_selected))]
  variant_calls_selected <- variant_calls_selected[, !grepl("->", 
                                                            colnames(variant_calls_selected))]
  out <- lapply(unique(patient), function(pa) {
    if (sum(mutation.bypatient == pa) == 0) 
      return(NULL)
    MN <- pullcounts.vars(BaseCounts[patient == pa], colnames(variant_calls_selected)[mutation.bypatient == pa])
    o <- mutationCallsFromMatrixCoverage(t(MN$M), t(MN$N))
  })
  names(out) <- unique(patient)
  out$blacklist <- rownames(varcount.bypatient[multipatient, ])
  out$blacklist <- gsub("(\\d+)(.+)", "\\1 \\2", out$blacklist)
  return(out)
}


# function to determine the CNA status of entire chromosome from the copykat output (trisomy or monosomy).
# returns a table with cell barcodes and the status of the input chromosomes (reference, trisomy and monosomy)
cna_full_chromosome <- function(copykat_output, chromosomes, trisomy = TRUE, valid_barcodes){
  
  
  status <- lapply(1:length(chromosomes), function(i){
    
    # average the value for all regions in the selected chromosome per cell
    table <- copykat_output$CNAmat %>% filter(chrom == chromosomes[i]) %>% 
                pivot_longer(cols = c(-chrompos, -abspos, -chrom),
                             names_to = "cell_barcode", values_to = "cna") %>% 
                group_by(cell_barcode) %>% 
                dplyr::summarise(cna = mean(cna)) %>% 
                mutate(cell_barcode = gsub("\\..+","", cell_barcode))
    
    # get cells with CNA
    if(trisomy[i]){
      
      table <- table %>% mutate(cna = ifelse(cna > 0.2, "trisomy", "reference"))
      
    }else{table <- table %>% mutate(cna = ifelse(cna < -0.2, "monosomy", "reference"))}
    
    # get dropout cells (copyKAT removes them)
    dropout_barcodes <- valid_barcodes[!valid_barcodes %in% table$cell_barcode]
    
    # add this as dropouts to the final table
    table <- table %>% bind_rows(data.frame(cna = "dropout", 
                                            cell_barcode = dropout_barcodes)) %>% 
              arrange(cell_barcode)
    
    # give the colnames correspongind to the selected chromosome
    colnames(table) <- c("cell_barcode", paste0("cnaChr", chromosomes[i]))
    
    # remove cell_barcode
    table %>% dplyr::select(-cell_barcode)
    
  })
  
  # bind columns for each CNA
  cna_merged <- do.call("cbind", status)
  
  
  # create a table with cell_barcodes and status for CNAs 
  final_table <- data.frame(cell_barcode = sort(valid_barcodes), cna_merged)
  
  return(final_table)
  
}


# function to add CNA status to mitoClone object
cna_mitoClone <- function(mitoClone_object, cna_table){
  
  
  for(i in length(colnames(cna_table))){
    
    # make M matrix by putting 10 reads to cells with CNA
    cna_m_matrix <- matrix(ifelse(cna_table[,colnames(cna_table)[i]] == "1", 10, 0),
                           dimnames = list(rownames(cna_table), colnames(cna_table)[i]))
   
    # make N matrix by putting 10 reads to cells with CNA 
    cna_n_matrix <- matrix(ifelse(cna_table[,colnames(cna_table)[i]] == "0", 10, 0),
                           dimnames = list(rownames(cna_table), colnames(cna_table)[i]))
    
    # create ternary matrix
    cna_ternary <- matrix(cna_table[,colnames(cna_table)[i]],
                          dimnames = list(rownames(cna_table), colnames(cna_table)[i]))
    
    # add entry so that the CNA is used for clustering
    cna_cluster <- setNames(TRUE, colnames(cna_table)[i])
    
    
    # add matrices to mitoClone object
    mitoClone_object@ternary <- cbind(mitoClone_object@ternary, as.matrix(cna_ternary))
    mitoClone_object@M <- cbind(mitoClone_object@M, as.matrix(cna_m_matrix))
    mitoClone_object@N <- cbind(mitoClone_object@N, as.matrix(cna_n_matrix))
    mitoClone_object@cluster <- c(mitoClone_object@cluster, cna_cluster)
    

  }
  
  return(mitoClone_object)
}
  
  
# modification of quickcluster function which includes patient name as title and removes the column names
quick_cluster_mod <- function (mutcalls, binarize = F, drop_empty = T, sample, ...){
  if (drop_empty) 
    mutcalls@ternary <- mutcalls@ternary[apply(mutcalls@ternary, 
                                               1, function(x) any(x == "1")), ]
  if (binarize) 
    converted <- t(apply(mutcalls@ternary, 2, function(x) ifelse(x == 
                                                                   "1", 1, ifelse(x == "0", -1, 0))))
  if (!binarize) 
    converted <- t(mutcalls@M/(mutcalls@M + mutcalls@N))
    converted[is.na(converted)] = 0
  clustered <- pheatmap::pheatmap(converted[mutcalls@cluster, ], main = sample, show_colnames = F,
                                  treeheight_col = 0, treeheight_row = 0,
                                  na_col = "gray", ...)
}


# function to filter blacklist variants with association tests to cell populations
# by default keep variants which are associated to any of the cell types with an adj. p_value < 0.01
filter_blacklist_variants <- function(clones_object, cell_types_table, p_value = 2,
                                      outpath = "plots"){

  
  # make ternary matrix for cell types present in the sample
  cell_type_ternary <- ternary_cell_type(cell_types = cell_types_table %>% 
                                                pull(cell_type) %>% unique(),
                                              table = cell_types_table) 
  
  
  
  # Association tests with all populations
  fisher_tests <- association_tests(mitoClone_object = clones_object,
                                         ternary_table = cell_type_ternary)
  
  
  
  # make qplots for populations in which at least on variant is significantly associated
  qplots_fisher <- make_qplots(fisher_tests, cell_types = TRUE)
  
  
  # save association plots in the corresponding folder
  png(file.path(outpath, "fisher_cell_types_blacklist_variants.png"),
      height = 15, width = 10, res = 200, units = "in")
  grid.arrange(grobs = qplots_fisher, ncol = 3)
  dev.off()
  
  
  # get variants by choosing those with an adjusted p_value < 0.01
  mito_variants <- do.call(rbind, fisher_tests) %>% 
    filter(pval > p_value) %>% 
    pull(name) %>% 
    unique()
  
  
  # get variants not associated with any population
  non_sig_variants <- colnames(clones_object@ternary)[which(!colnames(clones_object@ternary) %in% mito_variants)]
  
  
  # do not use these variants for clustering
  clones_object@cluster[non_sig_variants] <- FALSE
  
  
  # return mitoClone object
  return(clones_object)
  
  
}


# function to get WTA attributes of single cells from seurat object
get_wta_features <- function(sample_name, barcodes_table, 
                             seurat_path = "../../../gene_expression/data"){
  
  # make Seurat object for AKLA1b
  count_matrix_rna <- Read10X(data.dir = file.path(seurat_path,sample_name,
                                                   "outs/filtered_feature_bc_matrix/"))$`Gene Expression`
  
  # build the corresponding Seurat object. 
  seurat <- CreateSeuratObject(counts = count_matrix_rna, assay = "RNA")
  
  
  # filter out removed barcodes after doublet detection
  seurat <- seurat[,paste0(barcodes_table %>% filter(sample == sample_name) %>% 
                                           pull(cell_barcode), "-1")]
  
  # compute percentage of mito reads/cell
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
  
  
  
  # extract metadata table
  features <- seurat@meta.data %>% as.data.frame() %>% 
                rownames_to_column(var = "cell_barcode") %>% 
                dplyr::select(-orig.ident) %>% 
                mutate(sample = sample_name)
  
  return(features)
  
}


# function to append coverage and mutation status of mitochondrial variants to a wta table
get_mito_features <- function(wta_table, clones_object_cohort, clones_object_blacklist, 
                              sample_name,
                              count_table_path = "../../data"){
  
  
  # get mutations labelled as CLUSTER = TRUE in either the Cohort or Blacklist call
  mutations <- c(gsub("^(X\\d.+)([A-Z]{1}.+)","\\1\\.\\2",
                      names(clones_object_cohort@cluster)[which(clones_object_cohort@cluster)]),
                 names(clones_object_blacklist@cluster)[which(clones_object_blacklist@cluster)]) %>% 
                unique()
  
  # positions of the mutations of interest
  positions <- gsub("^X(\\d+).+","\\1",mutations) %>% as.integer()
  
  
  # load output of mGATK for AKLA1b
  mgatk <- readRDS(file.path(count_table_path, sample_name,"mito_counts/count_table/final", 
                             paste0(sample_name, ".rds")))
  
  
  # get coverage/cell for positions of interest 
  coverage_mutations <- mgatk@assays@data[["coverage"]][positions, ] %>% t()
  
  # put mutation names as column names and cell_barcodes as rownames
  colnames(coverage_mutations) <- mutations
  rownames(coverage_mutations) <- colnames(mgatk)
  
  
  # make exception for P266
  if(sample_name == "P266_TAPseq"){sample_name <- "P266"}
  
  
  # add coverage info to features table
  features <- wta_table %>% filter(sample == sample_name) %>% 
                left_join(coverage_mutations %>% as.data.frame() %>% 
                rownames_to_column(var = "cell_barcode") %>% 
                pivot_longer(cols = -cell_barcode,
                             values_to = "coverage", names_to = "mutation"))
  
  # add alelle frequencies (i add an arbitrary number to avoid dividing by 0)
  af_cohort <- clones_object_cohort@M/(clones_object_cohort@M + clones_object_cohort@N + 0.0000001)
  
  # the cohort mutations lack a '.' in their name...
  colnames(af_cohort) <- gsub("^(X\\d.+)([A-Z]{1}.+)","\\1\\.\\2",colnames(af_cohort))
  
  # get mutations unique in cohort
  mut_cohort <- gsub("^(X\\d.+)([A-Z]{1}.+)","\\1\\.\\2", 
                     names(clones_object_cohort@cluster)[which(clones_object_cohort@cluster)])

  # get allele frequency for blacklist
  af_blacklist <- clones_object_blacklist@M/(clones_object_blacklist@M + clones_object_blacklist@N + 0.0000001)
  
  # get mutations which are unique to the blacklist call
  mut_blacklist <- colnames(af_blacklist)[which(!colnames(af_blacklist) %in% colnames(af_cohort) &
                                                  colnames(af_blacklist) %in% mutations)]
  
  
  features_cohort <- features %>% filter(mutation %in% mut_cohort) %>% 
                left_join(af_cohort[, mut_cohort] %>% as.data.frame() %>% 
                            rownames_to_column(var = "cell_barcode") %>%
                            mutate(cell_barcode = gsub("_.+", "-1", cell_barcode)) %>% 
                            pivot_longer(cols = -cell_barcode,
                                         values_to = "af", names_to = "mutation")) %>% 
                left_join(clones_object_cohort@ternary %>% as.data.frame() %>% 
                            rownames_to_column(var = "cell_barcode") %>%
                            mutate(cell_barcode = gsub("_.+", "-1", cell_barcode)) %>% 
                            pivot_longer(cols = -cell_barcode,
                                         values_to = "mutation_status", names_to = "mutation") %>% 
                            mutate(mutation = gsub("^(X\\d.+)([A-Z]{1}.+)","\\1\\.\\2", mutation))) %>% 
                as_tibble() %>% 
                mutate(mutation_status = ifelse(mutation_status == "1", "mutant",
                                                ifelse(mutation_status == "0", "reference", 
                                                       "dropout")),
                       mutation = as.factor(mutation),
                       coverage = coverage + 0.000001)
              
    # add mutations unique of Blacklist call  
    if(length(mut_blacklist) > 0){
      
      
      # get allele frequency for blacklist
      af_blacklist <- clones_object_blacklist@M/(clones_object_blacklist@M + clones_object_blacklist@N + 0.0000001)
      
      af_blacklist <- as.matrix(af_blacklist[,mut_blacklist])
      colnames(af_blacklist) <- mut_blacklist

      
      features_blacklist <- features %>% filter(mutation %in% mut_blacklist) %>% 
                left_join(af_blacklist %>% as.data.frame() %>% 
                            dplyr::select(mut_blacklist) %>% 
                            rownames_to_column(var = "cell_barcode") %>%
                            mutate(cell_barcode = gsub("_.+", "-1", cell_barcode)) %>% 
                            pivot_longer(cols = -cell_barcode,
                                         values_to = "af", names_to = "mutation")) %>% 
                left_join(clones_object_blacklist@ternary %>% as.data.frame() %>% 
                            rownames_to_column(var = "cell_barcode") %>%
                            mutate(cell_barcode = gsub("_.+", "-1", cell_barcode)) %>% 
                            pivot_longer(cols = -cell_barcode,
                                         values_to = "mutation_status", names_to = "mutation")) %>% 
                as_tibble() %>% 
                mutate(mutation_status = ifelse(mutation_status == "1", "mutant",
                                                ifelse(mutation_status == "0", "reference", 
                                                       "dropout")),
                       mutation = as.factor(mutation),
                       coverage = coverage + 0.000001)
      
      features <- bind_rows(features_cohort, features_blacklist)
      
    }else{features <- features_cohort}
  
  return(features)
                
}


# add coverage threshold to mutationCalls function. The default will be 10 reads
mutationCallsFromMatrixAdj <- function (M, N, cluster = NULL, metadata = data.frame(row.names = rownames(M)), 
                                             binarize = 0.01, min_coverage = 10) 
{
  colnames(M) <- make.names(colnames(M))
  colnames(N) <- make.names(colnames(N))
  binfun <- function(M, N) {
    alleleRatio <- M/(M + N)
    coverage <- M+N
    ternary_list <- lapply(1:ncol(M), function(x){
      
      data <- data.frame(af = alleleRatio[,x], cov_site = coverage[,x]) %>% 
        mutate(ternary = ifelse(cov_site < min_coverage, "?", ifelse(af>binarize, "1", "0"))) %>% 
        dplyr::select(ternary)
    })
    
    ternary_matrix <- do.call("cbind", ternary_list)
    colnames(ternary_matrix) <- colnames(M)
    
    return(as.matrix(ternary_matrix))
    
  }
  out <- new("mutationCalls", M = M, N = N, metadata = metadata, 
             ternary = binfun(M, N))
  if (!is.null(cluster)) 
    out@cluster <- cluster
  else {
    out@cluster <- apply(out@ternary != "?", 2, mean) > 
      0.2
  }
  out
}


mutationCallsFromBlacklistAdj <- function (BaseCounts, lim.cov = 20, min.af = 0.2, min.num.samples = 0.01 * 
                                                  length(BaseCounts), min.af.universal = min.af, universal.var.cells = 0.95 * 
                                                  length(BaseCounts), blacklists.use = blacklists, max.var.na = 0.5, 
                                                max.cell.na = 0.95, ...) 
{
  varaf <- parallel::mclapply(BaseCounts, function(x) {
    x <- x[, 1:4]
    zeroes <- rowSums(x) < lim.cov
    x.af <- x/(x + apply(x, 1, max))
    x.af <- reshape2::melt(x.af)
    colnames(x.af) <- c("pos", "nt", "af")
    x.af <- x.af[!(mito.dna[x.af$pos] == x.af$nt), ]
    x.af <- x.af[!(mito.dna[x.af$pos] == "N"), ]
    x.af$name <- paste0(x.af$pos, " ", mito.dna[x.af$pos], 
                        ">", x.af$nt)
    x.af$af[x.af$pos %in% which(zeroes)] <- NA
    x <- x.af$af
    names(x) <- x.af$name
    return(x)
  }, mc.cores = 10)
  varaf <- do.call(cbind, varaf)
  varaf <- varaf[rowSums(varaf > min.af, na.rm = TRUE) >= 
                   min.num.samples, ]
  is.names <- sapply(blacklists.use, function(x) typeof(x) == 
                       "character")
  # if (sum(is.names) > 0) {
  #   removal.names.list <- unique(unlist(blacklists.use[is.names]))
  #   varaf <- varaf[!row.names(varaf) %in% removal.names.list, 
  #                  ]
  # }
  # if (sum(!is.names) > 0) {
  #   removal.ranges.list <- unique(unlist(GenomicRanges::GRangesList(blacklists.use[!is.names])))
  #   varaf <- varaf[-c(S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2gr(row.names(varaf)), 
  #                                                                      removal.ranges.list))), ]
  # }
  varaf <- varaf[rowSums(varaf, na.rm = T) > 0, ]
  varaf <- varaf[!rowSums(varaf >= min.af.universal, na.rm = TRUE) >= 
                   universal.var.cells, ]
  varaf <- varaf[rowSums(is.na(varaf)) < max.var.na * NCOL(varaf), 
                 ]
  varaf <- varaf[, colSums(is.na(varaf)) < max.cell.na * NROW(varaf)]
  MN <- pullcounts.vars(BaseCounts, rownames(varaf), colnames(varaf))
  mutationCallsFromMatrixAdj(t(MN$M), t(MN$N), ...)
}


# function to reformat UMI counts for mutations of interest from a pickle object to an R tibble
pickle_to_df <- function(pickle_path = "data/", 
                         sample_name){
  
  # load the object
  pickle <- pd$read_pickle(paste0(pickle_path, sample_name, "/mito_counts/umi_counts_mutations.pickle"))
  
  # for each mutated site create a dataframe with cell_barcode, nUMIs, nreads/UMI as columns
  lapply(1:length(pickle), function(i){
    
    # subset data for the mutation of interest
    data <- pickle[[i]]
    
    # get number of UMIs/cell
    table <- data %>% lengths() %>% as.matrix() %>% 
              as.data.frame() %>% 
              rownames_to_column(var = "cell_barcode")
    
    colnames(table) <- c("cell_barcode", "umis")
    
    
    # get number of reads/UMI and cell
    final_table <- table %>% left_join(data %>% unlist() %>% 
                  as.data.frame() %>% 
                  rownames_to_column(var = "cell_barcode") %>%
                  separate(cell_barcode, into = c("cell_barcode", "umi"), sep = "\\.") %>% 
                  group_by(cell_barcode) %>% 
                  dplyr::summarise(reads_umi = mean(.))) %>% as_tibble() 
    
  })
  
}

ternary_tcells <- function(cell_type_table,
                           tcells = c("CD4+ memory T cells", "CD8+CD103+ tissue resident memory T cells",
                                      "CD4+ naive T cells", "CD56dimCD16+ NK cells", "CD8+ effector memory T cells",
                                      "CD8+ naive T cells", "CD8+ central memory T cells", 
                                      "CD56brightCD16- NK cells", "GammaDelta T cells",
                                      "CD4+ cytotoxic T cells", "NK T cells", "CD69+PD-1+ memory CD4+ T cells")){
  
  
  
  new_table <- cell_type_table %>% 
    mutate(tcell = ifelse(celltype %in% tcells, "1", "0"),
           non_tcell = ifelse(celltype %in% tcells, "0", "1"))
  
  return(new_table)
  
}
  

# function to get mitochondrial M and N matrices for a selected subset of variants
mito_matrix <- function(seurat, variants, mito_path = "~/cluster/project/AML/mito_mutations/data/",
                        name = "", timepoints = FALSE, samples = "", runs = c("first", "second", "third")){
  
  if(timepoints){
    
    # iterate through timepoints
    time_list <- lapply(1:length(samples), function(i){
      
      samp <- samples[i]
      
      # load mito count table
      mito_object <- readRDS(paste0(mito_path, samp, "/mito_counts/count_table/", samp, "_count_table.rds" ))
      
      # get cellbarcodes from timepoint of interest
      cb <- intersect(seurat@meta.data %>% filter(run == runs[i]  & patient == str_split(samp, "")[[1]][4]) %>% rownames(),
                      names(mito_object))
      
      # select only cells present in Seurat object
      mito_object <- mito_object[cb]
      
      # iterate through the selected variants  
      var_list <- lapply(variants, function(j){
        
        # get position and reference and alternative alleles
        pos <- as.integer(gsub("X([0-9]+).+$", "\\1", j))
        ref <- gsub(".+([A-Z]{1}).+", "\\1", j)
        alt <- gsub(".+([A-Z])$", "\\1", j)
        
        # get ref and alt for each single cell
        cell_counts <- lapply(names(mito_object), function(k) mito_object[[k]][pos,c(ref,alt)])
        
        mat <- do.call("rbind", cell_counts)
        
        mat
        
        
      })
      
      # generate M and N matrices for the specific time point
      N_time <- do.call("cbind", lapply(var_list, function(i) i[,1]))
      M_time <- do.call("cbind", lapply(var_list, function(i) i[,2]))
      
      # add column and rownames
      rownames(N_time) <- rownames(M_time) <- names(mito_object)
      colnames(N_time) <- colnames(M_time) <- variants
      
      list(N = N_time, M = M_time)
      
    })
    
    # generate final M and N matrices
    object <- do.call(Map, c(f = rbind, time_list))
    
    return(object)
    
  }else{
    
    # get cellbarcodes from timepoint of interest
    cb <- colnames(seurat)
    
    # load mito count table
    mito_object <- readRDS(paste0(mito_path, name, "/mito_counts/count_table/", name, "_count_table.rds" ))
    
    # get cellbarcodes from timepoint of interest
    cb <- intersect(cb, names(mito_object))
    
    # select only cells present in Seurat object
    mito_object <- mito_object[cb]
    
    # iterate through the selected variants  
    var_list <- lapply(variants, function(j){
      
      # get position and reference and alternative alleles
      pos <- as.integer(gsub("X([0-9]+).+$", "\\1", j))
      ref <- gsub(".+([A-Z]{1}).+", "\\1", j)
      alt <- gsub(".+([A-Z])$", "\\1", j)
      
      # get ref and alt for each single cell
      cell_counts <- lapply(names(mito_object), function(k) mito_object[[k]][pos,c(ref,alt)])
      
      mat <- do.call("rbind", cell_counts)
      
      
    })
    
    # generate M and N matrices
    N <- do.call("cbind", lapply(var_list, function(i) i[,1]))
    M <- do.call("cbind", lapply(var_list, function(i) i[,2]))
    
    # add column and rownames
    rownames(N) <- rownames(M) <- names(mito_object)
    colnames(N) <- colnames(M) <- variants
    
    return(list(M = M, N = N))
    
    
  }
  
}
  