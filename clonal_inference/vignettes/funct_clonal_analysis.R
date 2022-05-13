# scripts with functions to carry out clonal analysis of patients
library(RColorBrewer)
library(patchwork)
library(viridis)
library(scmap)
library(SingleCellExperiment)
library(parallel)
library(aws.s3)

warriors_yellow <- "#FDBB30"
wine_red <- "#860038"
lakers_purple <- "#552583"
apple_gray <- "#7D7D7D"
jazz_blue <- "#002F6C"
bucks_green <- "#274E37"
dallas_blue <- "#2A6BAC"
bulls_red <- "#CE1141"
facebook_blue <- "#4267B2"
google_green <- "#0F9D58"
google_blue <- "#4285F4"
google_yellow <- "#F4B400"
google_red <- "#DB4437"
mitochondria_purple <- brewer.pal(5, "Set1")[4]
other_colour <- brewer.pal(5, "Set1")[5]
orange_set1 <- brewer.pal(8, "Dark2")[2]
golden <- brewer.pal(8, "Dark2")[7]
salmon <- brewer.pal(3, "Accent")[3]
dropout_gray <- "#d3d3d3"

n <- 200
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
discrete_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
discrete_colors <- discrete_colors[9:n]

# function to project cells to reference map from Sergio's paper. It only requires a seurat object.
project_reference <- function(seurat, neighbours = 5,
                              path = "~/cluster/project/AML/gene_expression/data/Sergio_figshare/Healthy.rds"){
  
  # Normalise data before projecting
  seurat <- seurat %>% NormalizeData()
  
  # get normalised counts
  norm_counts <- seurat@assays$RNA@data
  
  # load seurat object from Sergio's data
  Healthy <- readRDS(path)
  
  # only a handful of genes are used for the projection. I subset those genes for mapping
  data_projection <- norm_counts[intersect(rownames(norm_counts), rownames(Healthy)), ]
  
  
  # find nearest neighbour in the map 
  sce_Culture <- SingleCellExperiment(assays = list(normcounts =  as.matrix(data_projection)))
  logcounts(sce_Culture) <- normcounts(sce_Culture)
  rowData(sce_Culture)$feature_symbol <- rownames(sce_Culture)
  
  sce_All <- SingleCellExperiment(assays = list(normcounts = as.matrix(Healthy@assays$BOTH@data[rownames(data_projection),])))
  
  logcounts(sce_All) <- normcounts(sce_All)
  # use gene names as feature symbols
  rowData(sce_All)$feature_symbol <- rownames(sce_All)
  # remove features with duplicated names
  sce_All <- sce_All[!duplicated(rownames(sce_All)), ]
  
  sce_Culture<-setFeatures(sce_Culture,features =  rownames(sce_Culture))
  sce_Culture <- indexCell(sce_Culture)
  
  sce_All<-setFeatures(sce_All,features =  rownames(sce_All))
  sce_All <- indexCell(sce_All)
  
  
  Culture_Map <- scmapCell(
    projection = sce_Culture,
    index_list = list(
      sce_All = metadata(sce_All)$scmap_cell_index
    ),
    w = neighbours)
  
  # load pseudotime data
  s3load("Pseudotime_progenitors.rda", bucket = "triana-shiny",base_url="s3.embl.de",
         key="triana-shiny",secret="ammdD0s8b61RHUnTOtzC8u8uc2GfOKAI",
         use_https=F,region = "", check_region = F, verbose = T)
  pt <- pt[rownames(pt) %in% colnames(Healthy),]
  ct2short <- c("curve1"="cDC","curve2"="B cells","curve3"="Myelocytes","curve4"="Erythroid","curve5"="Megakaryocte")
  colnames(pt)<-ct2short[colnames(pt)]
  pt.all <- matrix(0, nrow = ncol(Healthy), ncol = ncol(pt), dimnames = list(colnames(Healthy), colnames(pt)))
  pt.all[rownames(pt),] <- pt
  pt.all[is.na(pt.all)] <- 0
  pt.all[pt.all==0] <- NA
  
  #calculate Coordinates
  Calc<-function(id,cult,pst){
    u <- cult[,id]
    xcoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[,1][u]
    ycoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[,2][u]
    x=median(xcoords)
    y=median(ycoords)
    x_mean=mean(xcoords)
    y_mean=mean(ycoords)
    meandist <- mean(sqrt((xcoords -x)^2 + (ycoords - y)^2))
    sdx <- sd(xcoords-x)
    sdy <- sd(ycoords-y)
    ct.t=table(Idents(Healthy)[u])
    ct.t<- ct.t[order(ct.t,decreasing = T)]
    Prop<-prop.table(ct.t)
    Prop<- Prop[order(Prop,decreasing = T)]
    nearest <- Healthy@assays$BOTH@data[rownames(sce_Culture),u]
    query <- normcounts(sce_Culture)[,id]
    #anglenn <- apply(nearest,2,function(x) angle(x, query))
    cornn<- apply(nearest,2,function(x) cor(x, query))
    pst.projected<-t(apply(pst, 2, function(xxx)  mean(na.omit(xxx[u]))))
    data.frame(row.names = id, x = x, y =y,x_mean=x_mean,y_mean=y_mean, ct = names(ct.t)[1],prop=Prop[1],meandist=meandist, sdx=sdx, sdy=sdy, cor = mean(cornn),pst.projected)
  }
  
  mapped <- mclapply(colnames(Culture_Map$sce_All[[1]]), Calc, cult = Culture_Map$sce_All[[1]],pst=pt.all, mc.cores = 10)
  
  mapped <- do.call(rbind,mapped) %>% rownames_to_column(var = "cell_barcode") %>% 
    dplyr::select(cell_barcode, x, y, ct, cor, Myelocytes)
  
  colnames(mapped) <- c("cell_barcode", "umapx", "umapy", "celltype", "score", "pseudo_myel")
  
  return(mapped)
}

# function to add simple celltypes to Seurat object
add_simple_ct_P80 <- function(seurat, table){
  
  ct2simple <- c(
    "Myelocytes" = "Monocytes",
    "Classical Monocytes" = "Monocytes",
    "Late promyelocytes" = "Early myeloid",
    "CD4+ memory T cells" = "T cells",
    "Early promyelocytes" = "Early myeloid",
    "Early erythroid progenitor" = "Erythroid",
    "Nonswitched memory B cells" = "B cells",
    "Erythro-myeloid progenitors" = "Erythroid",
    "Lymphomyeloid prog" = "Immature",
    "CD4+ naive T cells" = "T cells",
    "Late erythroid progenitor" = "Erythroid",
    "Class switched memory B cells" = "B cells",
    "Aberrant erythroid" = "Erythroid",
    "HSCs & MPPs" = "Immature",
    "CD8+CD103+ tissue resident memory T cells" = "T cells",
    "Mature naive B cells" = "B cells",
    "Conventional dendritic cell 1" = "Dendritic",
    "Eosinophil-basophil-mast cell progenitors" = "Other",
    "CD56dimCD16+ NK cells" = "NK cells",
    "CD8+ effector memory T cells" = "T cells",
    "CD8+ central memory T cells" = "T cells",
    "Plasma cells" = "B cells",
    "CD56brightCD16- NK cells" = "NK cells",
    "Pre-B cells" = "B cells",
    "GammaDelta T cells" = "T cells",
    "CD8+ naive T cells" = "T cells",
    "CD11c+ memory B cells" = "B cells",
    "Conventional dendritic cell 2" = "Dendritic",
    "Plasmacytoid dendritic cell progenitors" = "Dendritic",
    "Pro-B cells" = "B cells",
    "CD69+PD-1+ memory CD4+ T cells" = "T cells",
    "Non-classical monocytes" = "Monocytes",
    "Pre-pro-B cells" = "Immature",
    "Plasmacytoid dendritic cells" = "Dendritic",
    "CD4+ cytotoxic T cells" = "T cells",
    "Megakaryocyte progenitors" = "Other",
    "NK T cells" = "T cells",
    "Monocyte-like blasts" = "Monocytes",
    "Mesenchymal cells_1" = "Other")
  
  seurat@meta.data <- seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>% 
                        left_join(table %>% dplyr::select(cell_barcode, celltype)) %>% 
                        column_to_rownames(var = "cell_barcode")
  
  seurat$ct <- ifelse(seurat$celltype %in% names(ct2simple), ct2simple[seurat$celltype], "Other")
  
  return(seurat)
  
  
}

# function to add simple celltypes to Seurat object
add_simple_ct <- function(seurat, table){
  
  ct2simple <- c(
    "Myelocytes" = "Monocytes",
    "Classical Monocytes" = "Monocytes",
    "Late promyelocytes" = "Early myeloid",
    "CD4+ memory T cells" = "T cells",
    "Early promyelocytes" = "Early myeloid",
    "Early erythroid progenitor" = "Erythroid",
    "Nonswitched memory B cells" = "B cells",
    "Erythro-myeloid progenitors" = "Erythroid",
    "Lymphomyeloid prog" = "Immature",
    "CD4+ naive T cells" = "T cells",
    "Late erythroid progenitor" = "Erythroid",
    "Class switched memory B cells" = "B cells",
    "Aberrant erythroid" = "Erythroid",
    "HSCs & MPPs" = "Immature",
    "CD8+CD103+ tissue resident memory T cells" = "T cells",
    "Mature naive B cells" = "B cells",
    "Conventional dendritic cell 1" = "Dendritic",
    "Eosinophil-basophil-mast cell progenitors" = "Other",
    "CD56dimCD16+ NK cells" = "NK cells",
    "CD8+ effector memory T cells" = "T cells",
    "CD8+ central memory T cells" = "T cells",
    "Plasma cells" = "Plasma cells",
    "CD56brightCD16- NK cells" = "NK cells",
    "Pre-B cells" = "B cells",
    "GammaDelta T cells" = "T cells",
    "CD8+ naive T cells" = "T cells",
    "CD11c+ memory B cells" = "B cells",
    "Conventional dendritic cell 2" = "Dendritic",
    "Plasmacytoid dendritic cell progenitors" = "Dendritic",
    "Pro-B cells" = "B cells",
    "CD69+PD-1+ memory CD4+ T cells" = "T cells",
    "Non-classical monocytes" = "Monocytes",
    "Pre-pro-B cells" = "Immature",
    "Plasmacytoid dendritic cells" = "Dendritic",
    "CD4+ cytotoxic T cells" = "T cells",
    "Megakaryocyte progenitors" = "Other",
    "NK T cells" = "T cells",
    "Monocyte-like blasts" = "Monocytes",
    "Mesenchymal cells_1" = "Other")
  
  
  
  seurat$ct <- ifelse(seurat$celltype %in% names(ct2simple), ct2simple[seurat$celltype], "Other")
  
  return(seurat)
  
  
}

# Plot fraction and number of cells for each simple celltype
barplot_ct <- function(seurat){
  
  data <- as.data.frame(table(seurat$ct))
  colnames(data) <- c("celltype", "ncells")
  data$pct <- data$ncells/sum(data$ncells)
  data <- data %>%
            mutate(celltype = ifelse(celltype == "Early myeloid", "Early\nmyeloid", as.character(celltype)),
                  celltype = factor(celltype, levels = c("Immature", "Early\nmyeloid", "Monocytes", "Dendritic",
                                                          "T cells", "NK cells", "Erythroid", "B cells", "Other")))
  
  barplot <- ggplot(data, aes(x = celltype, y = pct*100, fill = celltype)) +
              geom_bar(stat = "identity") +
              theme_classic() +
              scale_y_continuous(expand = c(0,0))+
              scale_x_discrete(expand = c(0,1)) +
              geom_text(data = data %>% filter(pct > 0.05), aes(label = ncells),
                        position = position_stack(vjust = 0.5), color = "white") +
              scale_fill_viridis_d(option = "inferno", end = 0.8, direction = -1) +
              ylab("Percentage of cells") +
              theme(axis.title.x = element_blank(),
                    legend.position = "none")
  
  return(barplot)
  
}

# function to determine the proportion of high- and low-quality cells in the sample (nFeatures v.s. percent.mt)
qc_plots <- function(seurat, patient, min_feat = 1500, max_mt = 15){
  
  
  # make a table with nFeatures and percent.mt and determine whether cells passed the QC thresholds
  table_qc <- seurat@meta.data %>% 
                mutate(quality = factor(ifelse(percent.mt < max_mt & nFeature_RNA > min_feat, "good", "poor"),
                                        levels = c("poor", "good")))
  
  
  # make scatter plot
  scatter_feat_mt <- ggplot(data = table_qc, 
                            aes(x = nFeature_RNA, y = percent.mt, color = quality)) +
                            geom_point(size = 0.7) +
                            scale_color_manual(values = c(dropout_gray, google_green)) +
                            theme_classic() +
                            ggtitle(paste0("Features: ", min_feat, " - percent.mt: ", max_mt, "%"))+
                            guides(color = guide_legend(override.aes = list(size = 3.5))) +
                            theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5),
                                  legend.position = "none") 
  
  agg_table <- table_qc %>% group_by(quality) %>% 
                tally() %>% 
                dplyr::mutate(prop = n/sum(n),
                              sample = patient)
  
  # make stacked barplot
  stack_bp <- ggplot(agg_table, 
                        aes(x = sample, y = prop*100, fill = quality))+
                    geom_bar(stat = "identity") +
                    theme_classic() +
                    scale_y_continuous(expand = c(0,0))+
                    scale_x_discrete(expand = c(0,0.5)) +
                    geom_text(data = agg_table %>% filter(prop > 0.05), 
                              aes(label= n),
                              position = position_stack(vjust = 0.5), color = "white") +
                    scale_fill_manual(values = c(dropout_gray, google_green)) +
                    ylab("Percentage of cells")+
                    theme(legend.position = "none",
                          axis.title.x = element_blank())
  
  
  # count cells by ct and quality
  agg_table_ct <- table_qc %>% group_by(ct, quality) %>% 
                tally() %>% 
                dplyr::mutate(prop = n/sum(n))
  
  # make stacked barplot
  stack_bp_ct <- ggplot(agg_table_ct, 
                     aes(x = ct, y = prop*100, fill = quality))+
                     geom_bar(stat = "identity") +
                    theme_classic() +
                    scale_y_continuous(expand = c(0,0))+
                    scale_x_discrete(expand = c(0,0.5)) +
                     geom_text(data = agg_table_ct %>% filter(prop > 0.05), 
                               aes(label= n),
                                   position = position_stack(vjust = 0.5), color = "white") +
                    scale_fill_manual(values = c(dropout_gray, google_green)) +
                    ylab("Percentage of cells")+
                    theme(legend.position = "bottom",
                          axis.title.x = element_blank())
  
  # join all plots together
  return((scatter_feat_mt + stack_bp + plot_layout(widths = c(3,1)))/(stack_bp_ct)) 
    
}

# function to make UMAP with nucelar variants discretised
umap_nucl <- function(mut_data, seurat){
  
  # get number of mutations to create plots
  nmut <- mut_data %>% pull(symbol) %>% unique()
  
  # create table with UMAP coordinates and mutational status + CNAs
  umap_mut_table <- seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>% 
                      dplyr::select(cell_barcode, umapx, umapy) %>% 
                    left_join(mut_data %>% 
                    pivot_wider(id_cols = cell_barcode,
                                names_from = "symbol", 
                                values_from = "status") %>% 
                    mutate_at(vars(nmut), function(x) ifelse(is.na(x), "dropout", x)))
  
  # create one UMAP plot per mutation mutation
  umap_mut_list <- lapply(1:length(nmut), function(i){
    
    
    umap_plot <- ggplot(umap_mut_table,
                        aes_string(x = "umapx", y = "umapy", colour = nmut[i])) +
      geom_point(data = umap_mut_table %>% filter(.[[nmut[i]]] == "dropout"),
                 size = 0.75, color = "gray") +
      geom_point(data = umap_mut_table %>% filter(.[[nmut[i]]] != "dropout"),
                 size = 0.75)+
      ggtitle(nmut[i]) +
      scale_color_manual(values = setNames(c("gray", bulls_red, dallas_blue),
                                           c("dropout", "mutant", "reference"))) +
      theme_classic() +
      theme(axis.title = element_text(size = 12),
            axis.text = element_text(size = 10, color = "black"),
            legend.title = element_blank(),
            plot.title = element_text(size = 18, hjust = 0.5),
            legend.text = element_text(size = 15),
            legend.position = "bottom") +
      guides(colour = guide_legend(override.aes = list(size = 3.5)))
  })
  
  
  # put plots side by side 
  arranged <- wrap_plots(umap_mut_list, guides = "collect") & theme(legend.position = "bottom")
  
  return(list(plot = arranged, data = umap_mut_table))
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

# function to find mito variants
mutationCallsFromBlacklistAdj <- function (BaseCounts, lim.cov = 20, min.af = 0.2, min.num.samples = 0.01 * 
                                             length(BaseCounts), min.af.universal = min.af, universal.var.cells = 0.95 * 
                                             length(BaseCounts), blacklists.use = blacklists, max.var.na = 0.5, 
                                           max.cell.na = 0.95, cores = 8, ...) 
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
  }, mc.cores = cores)
  varaf <- do.call(cbind, varaf)
  varaf <- varaf[rowSums(varaf > min.af, na.rm = TRUE) >= 
                   min.num.samples, ]
  is.names <- sapply(blacklists.use, function(x) typeof(x) == 
                       "character")
  varaf <- varaf[rowSums(varaf, na.rm = T) > 0, ]
  varaf <- varaf[!rowSums(varaf >= min.af.universal, na.rm = TRUE) >= 
                   universal.var.cells, ]
  varaf <- varaf[rowSums(is.na(varaf)) < max.var.na * NCOL(varaf), 
                 ]
  varaf <- varaf[, colSums(is.na(varaf)) < max.cell.na * NROW(varaf)]
  MN <- pullcounts.vars(BaseCounts, rownames(varaf), colnames(varaf))
  mutationCallsFromMatrixAdj(t(MN$M), t(MN$N), ...)
}

# make ternary matrix for cell population (to be used in association tests)
ternary_cell_type <- function(seurat){
  
  data <- data.frame(tcell = ifelse(seurat$ct %in% c("T cells", "NK cells"), 1, 0),
                     myeloid = ifelse(seurat$ct %in% c("Immature", "Early myeloid", 
                                                       "Monocytes", "Dendritic"), 1, 0)) %>% 
            as.matrix()

  # add cell barcodes as rownames
  rownames(data) <- colnames(seurat)
  
  
  return(data)
  
  
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
make_qplots <- function(fisher_output, cell_types = FALSE, pval_thr = 2){
  
  
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
          data=na.omit(fisher_output[[sig_columns[i]]] %>% filter(pval > pval_thr))) + 
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
# function to get nuclear M and N matrices for a particular sample
nuclear_matrix <- function(seurat, variants, path = "~/cluster/project/AML/nuclear_mutations/data/",
                           name = "", timepoints = FALSE, samples = "", runs = c("first", "second", "third"),
                           barcodes = NULL){
  
  if(timepoints){
    
    # iterate through timepoints
    time_list <- lapply(1:length(samples), function(i){
      
      samp <- samples[i]
      
      # load mito count table
      nucl_object <- readRDS(paste0(path, samp, "/mutation_counts/count_table/", samp, "_count_table.rds"))
      
      # get cellbarcodes from timepoint of interest
      cb <- intersect(seurat@meta.data %>% filter(run == runs[i]  & patient == str_split(samp, "")[[1]][4]) %>% rownames(),
                      nucl_object$cell_barcode)
      
      # N matrix
      N <- nucl_object %>% filter(cell_barcode %in% cb & symbol %in% variants) %>%
        mutate(ref = ifelse(is.na(ref), 0, ref)) %>%
        pivot_wider(id_cols = c(cell_barcode), values_from = ref, names_from = symbol) %>% 
        column_to_rownames(var = "cell_barcode") %>% as.matrix()
      
      # M matrix
      M <- nucl_object %>% filter(cell_barcode %in% cb & symbol %in% variants) %>%
        mutate(alt = ifelse(is.na(alt), 0, alt)) %>%
        pivot_wider(id_cols = c(cell_barcode), values_from = alt, names_from = symbol) %>% 
        column_to_rownames(var = "cell_barcode") %>% as.matrix()
      
      list(M = M, N = N)
      
    })
    
    # generate final M and N matrices
    object <- do.call(Map, c(f = rbind, time_list))
    
    return(object)
    
  }else{
    
    # load mito count table
    nucl_object <- readRDS(paste0(path, name, "/mutation_counts/count_table/", name, "_count_table.rds"))
    
    if (name == "P80"){
      
      nucl_object <- readRDS(paste0(path, name, "/mutation_counts/count_table/", name, "_count_table.rds")) %>% 
                        filter(symbol != "RPS29") %>% 
                        bind_rows(readRDS(paste0(path, "P80_WT/mutation_counts/count_table/P80_WT_count_table.rds")) %>% 
                                    filter(symbol == "RPS29"))
      
    }
    
    # get cellbarcodes from timepoint of interest
    cb <- intersect(colnames(seurat),
                    nucl_object$cell_barcode)
    
    # N matrix
    N <- nucl_object %>% filter(cell_barcode %in% cb & symbol %in% variants) %>%
      mutate(ref = ifelse(is.na(ref), 0, ref)) %>%
      pivot_wider(id_cols = c(cell_barcode), values_from = ref, names_from = symbol) %>% 
      column_to_rownames(var = "cell_barcode") %>% as.matrix()
    
    # M matrix
    M <- nucl_object %>% filter(cell_barcode %in% cb & symbol %in% variants) %>%
      mutate(alt = ifelse(is.na(alt), 0, alt)) %>%
      pivot_wider(id_cols = c(cell_barcode), values_from = alt, names_from = symbol) %>% 
      column_to_rownames(var = "cell_barcode") %>% as.matrix()
    
    return(list(M = M, N = N))
    
  }
  
}


# function to get number of reads in the chromosome of interest
getCNVStatus <- function(s, M = NULL, N = NULL, chromosomes, mito = F,
                         partial = F, positions = NULL, host = "https://uswest.ensembl.org") {
  
  # this is needed when the uswest site is not available
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  
  # connect to ensembl database
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = host)
  # get all human genes with chromosomes names and start and end positions
  genesV2 = biomaRt::getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position", "end_position"), 
                           filters = "hgnc_symbol", values = rownames(s), mart = mart)
  genesV2 <- subset(genesV2, !grepl("GL|CHR_",chromosome_name))
  g2c <- genesV2$chromosome_name; names(g2c) <- genesV2$hgnc_symbol
  
  # get counts from Seurat and subset for gene names for which we have genomic info
  counts <- GetAssayData(s, assay="RNA", slot = "counts")
  counts <- counts[rownames(counts) %in% names(g2c),]
  
  row2c <- g2c[rownames(counts)]
  
  count.by.chrom <- sapply(unique(row2c), function(chr) {
    colSums(counts[row2c == chr,,drop=F])
  })
  
  # check if counts from a chromosome region should be retrieved
  if(length(which(partial == T))){
    
    partial_counts <- sapply(which(partial == T), function(i){
      
        chr <- gsub('^(\\d+).+',"\\1",chromosomes[i])
        start <- positions[[i]][1]
        end <- positions[[i]][2]
      
        # get genes in the region of interest
        part_genes <- genesV2 %>% filter(chromosome_name == chr,
                                         start_position >= start,
                                         end_position <= end) %>% 
                        pull(hgnc_symbol)
        
        # get counts the partial chromosomal region
        counts <- colSums(counts[part_genes,,drop=F])
      
    })
    
    # bind partial chromosome counts to total chromosome count matrix
    count.by.chrom <- cbind(count.by.chrom, partial_counts)
    
    # add column names
    colnames(count.by.chrom) <- c(as.character(1:22), "X", "Y", "MT",chromosomes[which(partial == T)])
    
  }
  
  if (!is.null(M) & !is.null(N)){
    usecells <- intersect(rownames(N),colnames(s))
    add.M <- count.by.chrom[usecells,chromosomes, drop=F]
    add.N <- replicate(length(chromosomes)+1, s$nCount_RNA[usecells])[,-1,drop=F]
    M <- M[usecells,,drop=F]
    N <- N[usecells,,drop=F]
    
    colnames(add.M) <- paste0("chr",chromosomes) -> colnames(add.N)
    list(M = cbind(M, add.M), N = cbind(N, add.N))
  } else {
    data.frame(count.by.chrom[,chromosomes,drop=F], total = s$nCount_RNA)
  }
  
}

# plot UMAP of mito variants from mito_matrix list object
umap_mito <- function(mito_list, seurat){
  
    var_list <- lapply(1:ncol(mito_list$M), function(i){
      
      data.frame(cell_barcode = rownames(mito_list$M),
                 alt = mito_list$M[,i],
                 ref = mito_list$N[,i],
                 symbol = colnames(mito_list$M)[i]) %>% 
      left_join(seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>% 
                dplyr::select(celltype, umapx, umapy, cell_barcode)) %>% 
      mutate(status = ifelse(alt > 0, "mutant", ifelse(ref > 0, "reference", "dropout")))
      
    })
    
    data <- do.call("bind_rows", var_list)
    
    umap <- ggplot(data, aes(x = umapx, y = umapy, colour = status)) +
                      geom_point(size = 0.75) +
                      scale_color_manual(values = setNames(c("gray", bulls_red, dallas_blue),
                                                           c("dropout", "mutant", "reference"))) +
                      theme_classic() +
                      theme(axis.title = element_text(size = 12),
                            axis.text = element_text(size = 10, color = "black"),
                            legend.title = element_blank(),
                            plot.title = element_text(size = 18, hjust = 0.5),
                            legend.text = element_text(size = 15),
                            legend.position = "bottom") +
                      guides(colour = guide_legend(override.aes = list(size = 3.5))) +
                      facet_wrap(vars(symbol))
  
    return(list(mito_table = data, plot = umap))
}

# heatmap showing the allele frequency of mutations as well the ct identity (T cell or non-T-cell)
heatmap_vaf <- function(seurat, patient, mito_matrix = NULL, nucl_table = NULL, mito_selected = "", nucl_variants = "",
                        nucl = T, mito = F, cnv = F, chr = "", mid_point = 0.5, cluster_muts = T, cb_order = ""){
  
    object <- NULL
    
    # get N and M matrices
    # get nuclear M and N matrices
    if(nucl){
      
      nucl_matrix <- nuclear_matrix(seurat, name = patient, nucl_variants, timepoints = F, samples = "", runs = "")
      
      if(!is.null(mito_matrix)){
        
        # intersect common barcodes
        cb <- intersect(rownames(mito_matrix$N), rownames(nucl_matrix$N))
        
        nucl_names <- colnames(nucl_matrix$M)
        
        nucl_N <- as.matrix(nucl_matrix$N[cb,])
        nucl_M <- as.matrix(nucl_matrix$M[cb,])
        colnames(nucl_M) <- colnames(nucl_N) <- nucl_names
        
        mito_names <- colnames(mito_matrix$N)
        
        mito_N <- as.matrix(mito_matrix$N[cb,])
        mito_M <- as.matrix(mito_matrix$M[cb,])
        colnames(mito_M) <- colnames(mito_N) <- mito_names
        
        object <- list(M = cbind(nucl_M, mito_M),
                       N = cbind(nucl_N, mito_N))
        
      # in case there is only nuclear mutations 
      }else{
        
        object <- nucl_matrix
        cb <- rownames(nucl_matrix$M)
      
      }
      
    }
    
    # merge SNVs matrices
    if(is_null(object)){
      
      object <- mito_matrix
      
      cb <- rownames(mito_matrix$M)
      
    }
    
    # create matrix with VAF to plot in the heatmap
    heatmap_matrix <- object$M/(object$N+object$M)
    
    # add CNV counts if specified
    if(cnv){
      
      chr_counts <- getCNVStatus(seurat, chromosomes = chr)
      
      counts <- lapply(chr, function(i){
        
        ratio <- chr_counts[cb,paste0("X", i)]/chr_counts[,"total"]
        ratio <- data.frame(ratio = ratio/max(ratio), row.names = rownames(chr_counts)) 
        colnames(ratio) <- paste0("chr",i)
        ratio
        
      })
      
      # merge cnv data
      merged_counts <- do.call("cbind", counts) %>% as.matrix()
      
      # merge cnv data with mito and snv
      heatmap_matrix <- cbind(heatmap_matrix, merged_counts)
      
    }
    
    # rename mito mutations and add them as column names of the matrix
    if(mito){
      
      mutnames <- sapply(colnames(object$M), function(i){
        
        ifelse(grepl("\\.", i), mito_rename(i), i)
        
      })
      
    }else{mutnames <- colnames(object$M)}
    
    # add cnv names if present
    if(cnv){mutnames <- c(mutnames, chr)}
    
    
    names(mutnames) <- NULL
    
    colnames(heatmap_matrix) <- mutnames
    rownames(heatmap_matrix) <- cb
    
    # remove rows that have all NAs 
    heatmap_matrix <- heatmap_matrix[rowSums(is.na(heatmap_matrix)) != ncol(heatmap_matrix),]
    
    # annotate T cells and other celltypes
    ct_table <- seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>%  
                  dplyr::select(cell_barcode, ct) %>% 
                  dplyr::mutate(ct = factor(ifelse(ct %in% c("T cells", "NK cells"), "T cells", "Other"), 
                                      levels = c("T cells", "Other"))) %>% 
                  filter(cell_barcode %in% rownames(heatmap_matrix)) %>% 
                  arrange(ct)
    
    # create object with the heat map column annotations
    col_ann <- HeatmapAnnotation("Celltype" = ct_table$ct,
                                 col = list("Celltype" = setNames(c(wine_red, dallas_blue),
                                                                  c("Other", "T cells"))),
                                 simple_anno_size = unit(1, "cm"),
                                 annotation_legend_param = list(Celltype = list(ncol = 1, title_position = "topcenter")))
    
    heatmap_matrix <- heatmap_matrix[ct_table$cell_barcode,] %>% t()
    
    # convert NAs to 0 and 0 to 0.0001 and adjust colours accordingly (allows row clustering)
    heatmap_matrix <- ifelse(is.na(heatmap_matrix), 0, ifelse(heatmap_matrix == 0, 0.0001, heatmap_matrix))

    # make heatmap in which the order of the rows is given by the clonal information
    heatmap <- Heatmap(heatmap_matrix,
                       top_annotation = col_ann,
                       heatmap_legend_param = list(title = "VAF", legend_direction = "vertical",
                                                   legend_width = unit(3, "cm"),
                                                   title_position = "topcenter"),
                       col = colorRamp2(c(0,0.0001,mid_point,1), c(dropout_gray, jazz_blue, warriors_yellow, google_red)),
                       show_heatmap_legend = T,
                       cluster_columns = cluster_muts,
                       cluster_rows = T,
                       show_column_dend = F,
                       show_row_dend = F,
                       show_row_names = T,
                       height = unit(1.625*nrow(heatmap_matrix), "cm"),
                       show_column_names = F,
                       use_raster = F,
                       na_col = dropout_gray,
                       column_title = "Cells",
                       column_title_side = "bottom")

    draw(heatmap, merge_legend = TRUE)
    #draw(heatmap)
    
    return(heatmap)
    
  
}

# function to plot the number of mutant and reference cells for mito and nuclear mutations for each 
# of the cell types
barplot_muts <- function(mito_table = NULL, nucl_table = NULL, seurat, mito = T, nucl = T){
  
  if(mito & nucl){
    
    data <- bind_rows(mito_table %>% dplyr::select(cell_barcode, symbol, status),
                      nucl_table %>% filter(cell_barcode %in% colnames(seurat)) %>% 
                        dplyr::select(cell_barcode, symbol, status))
    
  }else if(mito){
    
    data <- mito_table %>% dplyr::select(cell_barcode, symbol, status)
    
  }else{data <- nucl_table %>% filter(cell_barcode %in% colnames(seurat)) %>% 
                  dplyr::select(cell_barcode, symbol, status)}
  
  data <- data %>% 
            left_join(seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>% 
                        dplyr::select(cell_barcode, ct)) %>% 
            group_by(symbol, ct, status) %>% 
            tally() %>% 
            ungroup() %>% 
            group_by(symbol) %>% 
            mutate(pct = n/sum(n)*100)
  
  barplot <- ggplot(data, aes(x = ct, y = pct, fill = status)) +
                      geom_bar(stat = "identity") +
                      theme_classic() +
                      scale_y_continuous(expand = c(0,1))+
                      scale_x_discrete(expand = c(0,1)) +
                      geom_text(data = data %>% filter(pct > 5), aes(label = n),
                                position = position_stack(vjust = 0.5), color = "white") +
                      scale_fill_manual(values = setNames(c("gray", bulls_red, dallas_blue, bulls_red, bulls_red),
                                         c("dropout", "mutant", "reference", "trisomy", "deletion"))) +
                      ylab("Percentage of cells") +
                      theme(axis.title.x = element_blank(),
                            legend.position = "bottom",
                            legend.title = element_blank(),
                            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                      facet_wrap(vars(symbol))
  
  return(list(data = data, plot = barplot))
  
}

# function to save M and N matrix as json to load in python 
# I also add the celltype as additional column
mutcalls_to_json <- function(seurat = NULL, patient = "", 
                             nucl = F, nucl_variants = NULL,
                             mito = T, mito_variants = NULL,
                             exome = T, bulk_mito = F,
                             CNV = FALSE, type_cnv = NULL,
                             CNV_celltypes = FALSE, output, 
                             chromosomes = NULL,
                             partial_cnv = F,
                             positions = NULL,
                             host = "https://uswest.ensembl.org",
                             CNV_priors = "~/cluster/lvelten/Analysis/AML/MethodsDev/clonalAssignment/ourcohort.partial.prior.rds",
                             cnv_mapper = "~/cluster/lvelten/Analysis/AML/MethodsDev/clonalAssignment/ourcohort.partial.celltypemapper.RDS",
                             time = F, time_sampl = "",
                             runs = c("first", "second", "third")){
  
  object <- NULL
  
  # get N and M matrices
  if(mito){
    
    mito_matrix_object <- mito_matrix(seurat, name = patient, mito_variants, timepoints = time, samples = time_sampl, runs = runs)
    
  }else{mito_matrix_object <- NULL}
  
  # get nuclear M and N matrices
  if(nucl){
    
    if(patient == "P19003-11"){
      patient = "P19003_1"}
    
    nucl_matrix <- nuclear_matrix(seurat, name = patient, nucl_variants, timepoints = time, samples = time_sampl, runs = runs)
    
    if(!is.null(mito_matrix_object)){
      
      # intersect common barcodes
      cb <- intersect(rownames(mito_matrix_object$N), rownames(nucl_matrix$N))
      
      nucl_names <- colnames(nucl_matrix$M)
      
      nucl_N <- as.matrix(nucl_matrix$N[cb,])
      nucl_M <- as.matrix(nucl_matrix$M[cb,])
      colnames(nucl_M) <- colnames(nucl_N) <- nucl_names
      
      mito_names <- colnames(mito_matrix_object$N)
      
      mito_N <- as.matrix(mito_matrix_object$N[cb,])
      mito_M <- as.matrix(mito_matrix_object$M[cb,])
      colnames(mito_M) <- colnames(mito_N) <- mito_names
      
      object <- list(M = cbind(nucl_M, mito_M),
                     N = cbind(nucl_N, mito_N))
      
      # in case there is only 
      
    }else{object <- nucl_matrix}
    
  }
  
  # merge SNVs matrices
  if(is_null(object)){
    
    object <- mito_matrix_object
    
  }
  
  # if there are CNVs get the number of counts in the chromosome of interest and the total number of counts
  if(CNV){
    
    # get counts on the chromosome of interest and nCount_RNA
    object <- getCNVStatus(seurat, M = object$M, N = object$N, chromosomes = chromosomes, partial = partial_cnv,
                           positions = positions, host = host) 
    
  }
  
  # add column names as item in the list
  object$colnames <- colnames(object$M)
  object$mut_names <- colnames(object$M)
  
  # add UMAP coordinates
  meta <- seurat@meta.data %>% rownames_to_column(var = "bc") %>% 
    filter(bc %in% rownames(object$M)) %>% 
    column_to_rownames(var = "bc")
  meta <- meta[rownames(object$M),]
  
  object$umapx <- meta$umapx
  object$umapy <- meta$umapy
  object$cell_barcode <- rownames(object$M)
  object$class_assign <- ifelse(meta$ct == "T cells", 0, 1)
  object$class_names <- c("Tcells", "myeloid")
  #object$timepoint <- ifelse(meta$day == "d0", 0, ifelse(meta$day %in% c("d15", "d21"), 1, 2))
  
  # add celltype and CNV h priors by celltype
  if(CNV_celltypes){
    
    # read cnv priors for all chromosomes
    cnv_priors <- readRDS(CNV_priors)
    
    # read celltype mapper (we have different groups of celltypes based on similarity of chromosome counts)
    mapper <- readRDS(cnv_mapper)
  
    # map celltypes to chromosome count clusters
    object$celltype <- as.integer(as.factor(ifelse(meta$celltype %in% names(mapper), mapper[meta$celltype], NA)))
    object$celltype_names <- mapper %>% unique() %>% sort()
    
    # rename chromosomes
    chrom <- gsub("^(\\d+).+", "\\1", chromosomes)
    
    # get cnv priors for the chromosomes of interest
    priors <- lapply(1:length(chrom), function(x){
      
      i <- chrom[x]
      
      # get partial priors if indicated
      if(partial_cnv[x]){
        
        i <- paste0(i, ":", positions[[x]][1], "-",positions[[x]][2])}
        
      list(mean = as.double(cnv_priors[[i]][["par"]][["mct"]]), 
           sd = cnv_priors[[i]][["par"]][["sct"]]) 
    })
    
    names(priors) <- chrom
    
    object$cnv_priors <- priors
    
  }
  
  # add bulk exome (if available)
  if(exome & nucl){
    
    # load output table from exome pipeline
    exome_out <- read_csv(paste0("../../exome/results/", patient, "/annotated_variants.csv"))
    
    list_exome <- lapply(nucl_variants, function(i){
      
      counts <- exome_out %>% filter(symbol == i) %>% 
                dplyr::mutate(counts_tcells = round(af_tcells*depth_tcells),
                              counts_myeloid = round(af_tumor*depth_tumor_cells))
      
      list_counts <- list(m = list(tcells = c(counts$counts_tcells), myeloid = c(counts$counts_myeloid)),
                          n = list(tcells = c(counts$depth_tcells-counts$counts_tcells), 
                                   myeloid = c(counts$depth_tumor_cells-counts$counts_myeloid)))
      
    })
    
    # create a list of M and N counts for T-cells and non-tcells
    bulk_nucl_m <- list(unlist(map(map(list_exome, "m"), "tcells")), unlist(map(map(list_exome, "m"), "myeloid")))
    bulk_nucl_n <- list(unlist(map(map(list_exome, "n"), "tcells")), unlist(map(map(list_exome, "n"), "myeloid")))
    
  }else{
    
    bulk_nucl_m <- bulk_nucl_n <- list(rep(0, length(nucl_variants)))
    
  }
  
  # add bulk counts from mitochondrial mutations (if available)
  # add bulk exome (if available)
  if(bulk_mito & mito){
    
    # load output table from exome pipeline
    bulk_mito <- read_csv(paste0("../../bulk_atac/data/", patient, "/variants/variants_table.csv"))
    
    list_mito <- lapply(mito_variants, function(i){
      
      counts <- bulk_mito %>% filter(name_mut == i) %>% 
        dplyr::mutate(counts_tcells = round(af_tcells*depth_tcells),
                      counts_myeloid = round(af_tumor*depth_tumor))
      
      list_counts <- list(m = list(tcells = c(counts$counts_tcells), myeloid = c(counts$counts_myeloid)),
                          n = list(tcells = c(counts$depth_tcells-counts$counts_tcells), 
                                   myeloid = c(counts$depth_tumor-counts$counts_myeloid)))
      
    })
    
    # create a list of M and N counts for T-cells and non-tcells
    bulk_mito_m <- list(unlist(map(map(list_mito, "m"), "tcells")), unlist(map(map(list_mito, "m"), "myeloid")))
    bulk_mito_n <- list(unlist(map(map(list_mito, "n"), "tcells")), unlist(map(map(list_mito, "n"), "myeloid")))
    
  }else{
    
    bulk_mito_m <- bulk_mito_n <- list(rep(0, length(mito_variants)), rep(0, length(mito_variants)))
    
  }
  
  # we might discover mito variants for which we don't have bulk atac counts, therefore we should fill 0s 
  if(length(bulk_mito_m[[1]]) < length(mito_variants)){
    
    miss_mito <- length(mito_variants)- length(bulk_mito_m[[1]])
    
    bulk_mito_m <- Map(c, bulk_mito_m, rep(0, miss_mito))
    bulk_mito_n <- Map(c, bulk_mito_n, rep(0, miss_mito))
    
  }
  
  # join bulk counts from nuclear and mito mutations
  if(nucl & mito){
    
    bulk_M <- Map(c, bulk_nucl_m, bulk_mito_m)
    bulk_N <- Map(c, bulk_nucl_n, bulk_mito_n)
  
  }else if(nucl){
    
    bulk_M <- bulk_nucl_m
    bulk_N <- bulk_nucl_n
    
  }else{
    
    bulk_M <- bulk_mito_m
    bulk_N <- bulk_mito_n
  }
  
  object$bulk_M <- bulk_M
  object$bulk_N <- bulk_N
  
  # add R prior and fill the bulk data for CNVs with 0s
  if(!CNV){
    
    object$r_cnv <- rep(0, length(c(nucl_variants, mito_variants)))
    
  }else{
    
    object$r_cnv <- c(rep(0, length(c(nucl_variants, mito_variants))),
                      ifelse(type_cnv == "amp", 1.5, 0.5))
    
    object$bulk_M[[1]] <- c(object$bulk_M[[1]], rep(0, length(type_cnv)))
    object$bulk_M[[2]] <- c(object$bulk_M[[2]], rep(0, length(type_cnv)))
    object$bulk_N[[1]] <- c(object$bulk_N[[1]], rep(0, length(type_cnv)))
    object$bulk_N[[2]] <- c(object$bulk_N[[2]], rep(0, length(type_cnv)))
    
  }
  
  # add h_alpha and h_beta priors
  object$h_alpha <- c(2,1000,1)
  object$h_beta <- c(100,1000,1)
  
  # add vector with mutation types
  object$mut_type <- c(rep(1, length(nucl_variants)), rep(2, length(mito_variants)), rep(0, length(type_cnv)))
  
  # save matrices as json
  write_json(object, paste0(output, patient, ".json"))
  
  return(object)
  
}

# function to make a dataframe with all information about a particular sample (after running the probabilistic model)
tidy_pickle <- function(pickle, clone_names, tree = "0"){
  
  data <- data.frame(cell_barcode = pickle$cell_barcode,
                     umapx = pickle$umapx, 
                     umapy = pickle$umapy,
                     clone_n = apply(pickle$clonal_prob[[tree]], 1, which.max),
                     clonal_probability = apply(pickle$clonal_prob[[tree]], 1, max)) %>%  
    mutate(clone = factor(clone_names[clone_n], levels = clone_names),
           leukemia_prob = 1-pickle$clonal_prob[[tree]][,1], 
           status = ifelse(clone_n != 1, "cancer", "healthy")) %>% 
    arrange(clone_n)
  
  return(data)
  
}

# function to get colours according to clonal probabilities
getColors <- function(prob, f= function(x) x, f2 = f, rot = 0, get = "color") {
  base <- seq(0,1,length.out = length(prob)+1)[-1] + rot
  base[base > 1] <- base[base > 1] -1
  posy <- 0.9*sin(2*pi * base)
  posx <- 0.9*cos(2*pi * base)
  newposx <- sum(prob*posx)
  newposy <- sum(prob*posy)
  r <- sqrt(newposx^2 + newposy^2)
  angle <- atan(newposy / newposx)
  angle <- ifelse(newposx < 0 , angle + pi,
                  ifelse(newposx > 0 & newposy < 0, angle+2*pi, angle))
  angle <- angle / (2*pi)
  angle[is.na(angle)] <- 0
  #if (angle < 0) angle <- angle + 1
  colors <- hsv(base, 1,1)
  out <- hsv(angle, f(r),f(r))
  if (get == "color") return(out)
  if (get == "hue") return(angle)
}

# function to fix mito names
mito_rename <- function(string){
  
  pos <- gsub("X([0-9]+).+$", "\\1", string)
  
  ref <- gsub(".+\\.([A-Z]).+", "\\1", string)
  
  alt <- gsub(".+\\.([A-Z])$", "\\1", string)
  
  name <- paste0("mt:", pos, ref, ">", alt)
  
}

# function to get 
get_tree_attr <- function(pickle, tree){
  
  tree <- as.integer(tree)
  
  # get index position of the tree
  tree_index <- pickle$tree_indices[tree]+1
  
  # get binary matrix and parents and children lists
  tree_mat <- pickle$trees[[tree]]
  colnames(tree_mat) <- pickle$mutations_tree
  children <- pickle$children[[tree_index]] 
  parents <- pickle$parents[[tree_index]]
  
  # change name of mito mutations
  colnames(tree_mat) <- ifelse(grepl("\\.", colnames(tree_mat)), mito_rename(colnames(tree_mat)), colnames(tree_mat))
  
  # vector for direction of edges and node names
  node_names <- node_order <- rep("", nrow(tree_mat))
  start_nodes <- end_nodes <- rep("", nrow(tree_mat)-1)
  
  for(node in 1:nrow(tree_mat)){
    
    # 1st node is always healthy clone
    if(node == 1){
      
      node_names[1] <- start_nodes[1] <- "Healthy"
      node_order[1] <- 1
      
      # 2nd node always hangs from the root
    }else if (node == 2){
      
      # get node immediately after root
      nxt_node <- order(lengths(parents))[2]
      
      # get mutations in the node
      node_name <- paste(colnames(tree_mat)[which(tree_mat[nxt_node,]-tree_mat[1,] == 1)], collapse = "\n")
      end_nodes[1] <- node_names[2] <- node_name
      node_order[2] <- nxt_node
      
      # search which node comes downstream      
    }else{
      
      # get next node in the tree
      nxt_node <- order(lengths(parents))[node]
      
      # get new mutations in the node
      node_name <- paste(colnames(tree_mat)[which(tree_mat[nxt_node,]-tree_mat[as.integer(node_order[node-1]),] == 1)], 
                         collapse = "\n")
      
      # add name of node and the index
      node_names[node] <- node_name
      node_order[node] <- nxt_node
      
      # get parent node
      node_parent <- tail(node_order[node_order %in% as.character(parents[[nxt_node]]+1)], n = 1)
      start_nodes[node-1]  <- node_names[which(node_order == node_parent)] 
      
      end_nodes[node-1] <- node_name
      
    }
  }
  
  tree_mat <- tree_mat[as.integer(node_order),]
  rownames(tree_mat) <- node_names
  
  return(list("start_nodes" = start_nodes, "end_nodes" = end_nodes,
              "tree_mat" = tree_mat, "clone_names" = node_names[order(node_order)]))
  
}

# function to plot all clonal hierarchies selected by the model
plot_trees <- function(pickle, tree = "all", clone_cols = T){
  
  if(tree == "all"){
    
    # colors in case clone colors are not to be used
    tree_cols <- discrete_colors
    
    # generate tree for all selected trees
    tree_list <- lapply(1:length(pickle$trees), function(i){
      
      # get node order and attributes to build the tree
      tree_att <- get_tree_attr(pickle, i)
      
      if(clone_cols){
        
        combs <- combn(nrow(tree_att$tree_mat), 1, function(x)replace(numeric(nrow(tree_att$tree_mat)),x,1))
        colors <- rev(apply(combs, 1, getColors))
        
      }else{colors <- tree_cols[i]}
      
      # specify height of nodes depending on the length of mutation names
      node_heights <- ifelse(grepl("\n", rownames(tree_att$tree_mat)), 
                             0.3*lengths(str_split(rownames(tree_att$tree_mat),pattern = "\n")), 0.5)
      
      # specify node width
      node_widths <- ifelse(nchar(rownames(tree_att$tree_mat)) > 11 & grepl("\n", rownames(tree_att$tree_mat)) == F,
                            nchar(rownames(tree_att$tree_mat))*0.1, 1.1)
      
      # create the tree
      # create nodes
      nodes <- create_node_df(
        n = nrow(tree_att$tree_mat),
        label = rownames(tree_att$tree_mat),
        style = "filled", 
        shape = "oval", 
        fixed_size = TRUE,
        fontname = "Helvetica",
        fontsize = 11, 
        fontcolor = '#FFFFFF',
        height = node_heights, 
        width = node_widths,
        color = colors, 
        fillcolor = colors)
      
      # get ids for each label
      ids <- setNames(nodes$id,nodes$label)
      
      # create edges
      edges <- create_edge_df(
        from = ids[tree_att$start_nodes],
        to = ids[tree_att$end_nodes],
        color = '#000000',
        rel = "leading_to")
      
      # create graph
      graph <- create_graph(
        nodes_df = nodes, 
        edges_df = edges,
        graph_name = paste0("tree",pickle$tree_indices[i]),
        attr_theme = "tb")
      
      out <- graph
      
      #out <- list("tree" = graph, "clone_names" = tree_att$clone_names,
      #            "tree_mat" = tree_att$tree_mat)
      
    })
    
    # create a graph series
    grph_series <- create_graph_series()
    for (t in 1:length(tree_list)){
      
      grph_series <- add_graph_to_graph_series(grph_series, tree_list[[t]])
      
    }

    return(grph_series)
    
  }else{
    
    # get node order and attributes to build the tree
    tree_att <- get_tree_attr(pickle, tree)
    
    if(clone_cols){
      
      combs <- combn(nrow(tree_att$tree_mat), 1, function(x)replace(numeric(nrow(tree_att$tree_mat)),x,1))
      colors <- rev(apply(combs, 1, getColors))
      
    }else{colors <- "black"}

    # specify height of nodes depending on the length of mutation names
    node_heights <- ifelse(grepl("\n", rownames(tree_att$tree_mat)), 
                     0.3*lengths(str_split(rownames(tree_att$tree_mat),pattern = "\n")), 0.5)
    
    # specify node width
    node_widths <- ifelse(nchar(rownames(tree_att$tree_mat)) > 11 & grepl("\n", rownames(tree_att$tree_mat)) == F,
                          nchar(rownames(tree_att$tree_mat))*0.1, 1.1)
    
    # create the tree
    # create nodes
    nodes <- create_node_df(
      n = nrow(tree_att$tree_mat),
      label = rownames(tree_att$tree_mat),
      style = "filled", 
      shape = "oval", 
      fixed_size = TRUE,
      fontname = "Helvetica",
      fontsize = 11, 
      fontcolor = "#FFFFFF",
      height = node_heights, 
      width = node_widths,
      color = colors, 
      fillcolor = colors)
    
    # get ids for each label
    ids <- setNames(nodes$id,nodes$label)
    
    # create edges
    edges <- create_edge_df(
                from = ids[tree_att$start_nodes],
                to = ids[tree_att$end_nodes],
                color = "#000000",
                rel = "leading_to")
    
    # create graph
    graph <- create_graph(
                nodes_df = nodes, 
                edges_df = edges,
                attr_theme = "tb")
    
    render_graph(graph)
    
    # return tree and tree attributes
    return(list("tree" = graph, "clone_names" = tree_att$clone_names,
                "tree_mat" = tree_att$tree_mat))
  }
}

# function to export trees into a merged pdf file
export_trees <- function(grph_series, outdir = "plots/", h = 600, w = 300){
  
  # create an individual pdf for each tree
  for (i in 1:count_graphs_in_graph_series(grph_series)){
    
    export_graph(grph_series$graphs[[i]], 
                 title = grph_series$graphs[[i]]$graph_info$graph_name, 
                 file_name = paste0(outdir,"/",i,"_tree.pdf"), 
                 file_type = "pdf",
                 height = h, width = w)
    
  }
  
  # get created pdf files
  pdfs <- list.files(path = outdir, pattern = "tree.pdf$", full.names = T)
  
  # create merged pdf file
  out <- pdf_combine(pdfs, output = paste0(outdir, "/trees_clonetracer.pdf"))
  
  # remove individual pdfs
  system(paste0("rm ", paste(pdfs, collapse = " ")))
  
}

# function to make line plot showing the ELBO
print_elbo <- function(pickle, trees = "all", first_iter = 100, exclude = "", max_elbo = 1e6){
  
  # select all trees
  if(trees == "all"){trees = as.character(1:length(pickle$potential_trees)-1)}
  
  # make dataframe with ELBO values over iterations (by default the 1st 100 iterations are not shown)
  elbo_data <- do.call("cbind",pickle$ELBO) %>% as.data.frame() %>% 
                mutate(iteration = 1:length(pickle$ELBO[[1]])) %>% 
                pivot_longer(cols = -iteration, names_to = "tree", 
                             values_to = "ELBO") %>% 
                filter(iteration > first_iter & tree %in% trees) %>% 
                filter(! tree %in% exclude) %>% 
                mutate(tree = paste0("tree ", tree)) %>% 
                filter(ELBO < max_elbo)
  
  # make line plot
  elbo_plot <- ggplot(elbo_data, aes(x = iteration, y = ELBO, color = tree))+
                  geom_line() + 
                  theme_classic() +
                  scale_color_manual(values = discrete_colors)
  
  return(elbo_plot)
  
  
}

# function to make heatmap using complex heatmap package
make_heatmap <- function(pickle, seurat, pat, tree, cols = "",
                         cnv_type = NULL, cnv_pos = NULL, cnv_matrix = NULL, middle_point = 0.5,
                         cvg = F, cvg_max = 10, cvg_mid = 5, clust_rows = F){
  

  
  # get order of nodes in the tree and clone names
  tree_att <- get_tree_attr(pickle, tree)
  
  # get clonal probs and name them
  probs <- pickle$clonal_prob[[tree]]
  colnames(probs) <- tree_att$clone_names
  probs <- probs[,rownames(tree_att$tree_mat)]
  rownames(probs) <- pickle$cell_barcode
  pickle$clonal_prob[[tree]] <- probs
  
  # clone names in the order of clonal probability matrix
  clone_names <- tree_att$clone_names
  clones_order <- rownames(tree_att$tree_mat)
  
  # put data into a dataframe for plotting and order by leukemia probability
  umap_table <- tidy_pickle(pickle,
                            clone_names = clones_order,
                            tree = tree) %>% 
                    mutate(patient = pat) %>% 
                    left_join(seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>%  
                                dplyr::select(cell_barcode, ct)) %>% 
                    dplyr::mutate(ct = factor(ifelse(ct %in% c("T cells", "NK cells"), "T cells", "Other"), 
                                              levels = c("T cells", "Other")),
                                  row = as.factor(dplyr::row_number()),
                                  order = ifelse(clone == "Healthy", clonal_probability, 1-clonal_probability),
                                  status = factor(status, levels = c("healthy", "cancer"))) %>% 
                    group_by(clone) %>%
                    dplyr::arrange(desc(clonal_probability), .by_group = TRUE) 
  
  probs <- probs[, rev(1:length(clone_names))] 
  
  # get colours
  if (length(cols) == 1){
    
    combs <- combn(length(clone_names), 1, function(x)replace(numeric(length(clone_names)),x,1))
    clone_cols_discr <- apply(combs, 1, getColors)
    
  }else{clone_cols_discr <- cols}
  
  # make barplot of clonal probabilities
  barplot <- anno_barplot(probs[umap_table$cell_barcode,],
                          gp = gpar(fill = clone_cols_discr, col = NA), 
                          bar_width = 1, height = unit(2, "cm"), which = "column", 
                          border = T, axis_param = list(at = c(0,0.25,0.5,0.75,1)))
  
  # get single-cell colours based on clonal probabilities
  clone_cols = apply(probs[umap_table$cell_barcode,], 1, getColors, rot = 0)
  
  # create object with the heat map column annotations
  col_ann <- HeatmapAnnotation("Clonal\ncall" = 1:nrow(probs),
                               "Clonal\nprobability" = barplot,
                               "Celltype" = umap_table$ct,
                               col = list("Clonal\ncall" = setNames(clone_cols, 1:nrow(probs)),
                                          "Celltype" = setNames(c(wine_red, dallas_blue),
                                                                c("Other", "T cells"))),
                               show_legend = c(F, F, T),
                               simple_anno_size = unit(1, "cm"),
                               annotation_legend_param = list(Celltype = list(ncol = 1, title_position = "topcenter")))
  
  # create matrix with VAF to plot in the heatmap
  heatmap_matrix <- pickle$M/(pickle$N+pickle$M)
  
  mutnames <- sapply(pickle$mutations_matrix, function(i){
    
    ifelse(grepl("\\.", i), mito_rename(i), ifelse(grepl("chr", i), paste0(cnv_type, "\n",i), i))
    
  })
  
  names(mutnames) <- NULL
  
  colnames(heatmap_matrix) <- mutnames
  rownames(heatmap_matrix) <- pickle$cell_barcode
  
  # Unit scale trisomy reads
  if (length(cnv_type) > 0){
    
    for(i in 1:length(cnv_type)){
      
      max <- max(pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]])
      
      scaled_ratio <- pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]]/max
      
      heatmap_matrix[,cnv_pos[i]] <- scaled_ratio
      
    }
  }
  
  heatmap_matrix <- heatmap_matrix[umap_table$cell_barcode,] %>% t()
  
  # make heatmap in which the order of the rows is given by the clonal information
  heatmap <- Heatmap(heatmap_matrix,
                     top_annotation = col_ann, 
                     heatmap_legend_param = list(title = "VAF", legend_direction = "vertical",
                                                 legend_width = unit(3, "cm"),
                                                 title_position = "topcenter"),
                     col = colorRamp2(c(0,middle_point,1), c(jazz_blue, warriors_yellow, google_red)),
                     show_heatmap_legend = T,
                     cluster_columns = F,
                     cluster_rows = clust_rows,
                     show_column_dend = F,
                     show_row_dend = F, 
                     show_row_names = T,
                     height = unit(1.625*nrow(heatmap_matrix), "cm"),
                     show_column_names = F,
                     use_raster = F,
                     na_col = dropout_gray,
                     column_title = "Cells", 
                     column_title_side = "bottom") 
  
  # add legend
  lgd_list <- list(Legend(labels = gsub("\n","-",rownames(tree_att$tree_mat)), 
                          legend_gp = gpar(fill = rev(clone_cols_discr)),
                          title = "Clone"))
  
  # plot heatmap showing total coverage underneath
  if(cvg){
    
    # create matrix with VAF to plot in the heatmap
    heatmap_matrix <- pickle$N+pickle$M
    
    mutnames <- sapply(pickle$mutations_matrix, function(i){
      
      ifelse(grepl("\\.", i), mito_rename(i), ifelse(grepl("chr", i), paste0(cnv_type, "\n",i), i))
      
    })
    
    names(mutnames) <- NULL
    
    colnames(heatmap_matrix) <- mutnames
    rownames(heatmap_matrix) <- pickle$cell_barcode
    
    heatmap_matrix <- heatmap_matrix[umap_table$cell_barcode,] %>% t()
    
    # put a cap on the maximum number of reads
    heatmap_matrix <- ifelse(heatmap_matrix > cvg_max, cvg_max, heatmap_matrix)
    
    if(clust_rows){
      
      row_order <- rownames(heatmap@matrix)[row_order(heatmap)]
      
    }else{row_order <- rownames(heatmap@matrix)}
    
    # make heatmap in which the order of the rows is given by the clonal information
    heatmap_cvg <- Heatmap(heatmap_matrix,
                           #top_annotation = col_ann, 
                           heatmap_legend_param = list(title = "Total coverage", legend_direction = "vertical",
                                                       legend_width = unit(3, "cm"),
                                                       title_position = "topcenter"),
                           col = colorRamp2(c(0,max(as.numeric(unlist(heatmap_matrix)))), 
                                            c("white", bucks_green)),
                           show_heatmap_legend = T,
                           cluster_columns = F,
                           row_order = row_order,
                           cluster_rows = F,
                           show_column_dend = F,
                           show_row_dend = F, 
                           show_row_names = T,
                           height = unit(1.625*nrow(heatmap_matrix), "cm"),
                           show_column_names = F,
                           use_raster = F,
                           na_col = dropout_gray,
                           column_title = "Cells", 
                           column_title_side = "bottom") 
    
    ht_list = heatmap %v% heatmap_cvg
    
    draw(ht_list, annotation_legend_list = lgd_list, merge_legend = TRUE)
    
    return(list(plot = heatmap, legend = lgd_list))
    
  }else{
    
    draw(heatmap, annotation_legend_list = lgd_list, merge_legend = TRUE)
    
    return(list(plot = heatmap, legend = lgd_list))
    
  }
}

# function to make UMAP with clones
umap_clones <- function(pickle, seurat, tree, name){
  
  # load Healthy cells from Sergio's paper
  Healthy <- readRDS("~/cluster/project/AML/gene_expression/data/Sergio_figshare/Healthy.rds")
  
  bckgr <- data.frame(x = Embeddings(Healthy, reduction = "MOFAUMAP")[,1],
                      y = Embeddings(Healthy, reduction = "MOFAUMAP")[,2])
  
  rm(Healthy)
  
  # get order of nodes in the tree and clone names
  tree_att <- get_tree_attr(pickle, tree)
  
  clone_names <- tree_att$clone_names
  
  # get clonal probs and name them
  probs <- pickle$clonal_prob[[tree]]
  colnames(probs) <- tree_att$clone_names
  probs <- probs[,rownames(tree_att$tree_mat)]
  rownames(probs) <- pickle$cell_barcode
  probs <- probs[, rev(1:length(clone_names))] 
  pickle$clonal_prob[[tree]] <- probs
  
  # clone names in the order of clonal probability matrix
  clones_order <- rownames(tree_att$tree_mat)
  
  # get colors according to clonal probs
  clonal_cols <- data.frame(cell_barcode = pickle$cell_barcode,
                            color = apply(probs, 1, getColors),
                            hue = apply(probs, 1, getColors, get="hue"))
  
  data <- tidy_pickle(pickle,
                      clone_names = clones_order,
                      tree = tree) %>% 
                    dplyr::mutate(patient = name) %>% 
                    left_join(seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>%  
                                dplyr::select(cell_barcode, celltype, ct)) %>% 
                    dplyr::mutate(ct = factor(ct, levels = rev(c("Immature", "Early myeloid", "Erythroid", "Monocytes", "Dendritic",
                                                                 "Plasma cells","B cells", "NK cells", "T cells", "Other")))) %>% 
                    left_join(clonal_cols) %>% 
                    group_by(clone) %>%
                    dplyr::arrange(desc(clonal_probability), .by_group = TRUE) %>% 
                    ungroup() %>% 
                    dplyr::mutate(heat_order = row_number()) %>%
                    ungroup() 
  
  # Make heatmap
  # order cells by leukemia probability
  data <- data %>%
    group_by(clone) %>%
    dplyr::arrange(desc(clonal_probability), .by_group = TRUE) 
  
  # get discrete colors of clones
  combs <- combn(length(names), 1, function(x)replace(numeric(length(names)),x,1))
  clone_cols_discr <- rev(apply(combs, 1, getColors))
  
  
  # get single-cell colours based on clonal probabilities
  clone_cols = data.frame(cell_barcode = data$cell_barcode,
                          clone_cols = apply(probs[data$cell_barcode,], 1, getColors, rot = 0)) %>% 
    left_join(data %>% dplyr::select(cell_barcode, ct, clone, leukemia_prob)) %>% 
    dplyr::mutate(ct_simple = case_when(ct == "Immature" ~ "Immature\nmyeloid",
                                        ct %in% c("Monocytes", "Dendritic") ~ "Mature\nmyeloid",
                                        ct == "Early myeloid" ~ "Early\nmyeloid",
                                        ct == "Erythroid" ~ "Erythroid",
                                        ct %in% c("T cells", "NK cells") ~ "T & NK cells",
                                        ct %in% c("B cells", "Plasma cells") ~ "B cells", 
                                        ct == "Other" ~ "Other"),
                  value = row_number()) %>%
    group_by(ct_simple) %>% 
    dplyr::mutate(order = row_number(), 
                  total = 1) %>% 
    dplyr::mutate(ct_simple = factor(ct_simple, levels = rev(c("Immature\nmyeloid", "Early\nmyeloid",
                                                               "Mature\nmyeloid", "T & NK cells",  
                                                               "Erythroid", "B cells", "Plasma cells","Other")))) 
  
  
  # stacked barplot filled by clonal call
  clonal_mixture <- ggplot(clone_cols %>% ungroup(),
                           aes(x = ct_simple, y = total, fill = clone_cols)) + 
                          geom_bar(stat = "identity") + 
                          scale_fill_identity() +
                          scale_x_discrete(expand = c(0,0,0.01,0), position = "top") +
                          scale_y_continuous(expand = c(0,0,0.01,0)) +  
                          theme_classic() +
                          coord_flip()+
                          ylab("Number of cells") +
                          theme(legend.position = "none",
                                axis.title.y = element_blank(),
                                axis.text.x = element_text(size = 8.5),
                                axis.text.y = element_text(size = 11, color = "black"),
                                panel.border = element_rect(color = "black", fill = NA), 
                                strip.background = element_blank(), 
                                strip.text = element_blank()) 
  
  # UMAP with reference as background
  umap <- ggplot(data %>% arrange(cell_barcode),
                 aes(x = umapx, y = umapy, color = color)) +
                  geom_point(data = bckgr, aes(x = x, y = y, color = "#D5D5D5"), size = 0.1) +
                  geom_point(size = 0.45) +
                  theme_classic() +
                  scale_color_identity() +
                  coord_cartesian(xlim = c(-13, 15), ylim = c(-12, 12.5)) +
                  theme(axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.line = element_blank(),
                        panel.border = element_rect(color = "black", fill = NA),
                        axis.ticks = element_blank(),
                        legend.title = element_blank(),
                        legend.position = "none",
                        plot.margin=grid::unit(c(0.75,0,0,0), "mm"),
                        strip.text = element_text(size = 13))
  
  # combine plots
  umap/clonal_mixture + plot_layout(heights = c(1.75,1))
  
  # get single-cell colours based on clonal probabilities
  clone_cols = apply(probs[data$cell_barcode,], 1, getColors, rot = 0)
  
  # heatmap with clonal probabilities mixture colors
  umap_legend <- ggplot(data.frame(x = 1:nrow(probs), y = "Clonal\nmixture", color = clone_cols),
                        aes(x = x, y = y, fill = color)) + 
                        geom_tile() + 
                        scale_fill_identity() +
                        scale_x_discrete(expand = c(0,0)) +
                        scale_y_discrete(expand = c(0,0), position = "right", label = "Clonal call") +  
                        theme_classic() +
                        theme(axis.title = element_blank(),
                              axis.text = element_blank(),
                              axis.line = element_blank(),
                              legend.position = "none",
                              panel.border = element_rect(color = "black", fill = NA),
                              #axis.text.y = element_text(color = "black", size = 20),
                              axis.ticks = element_blank(),
                              plot.margin=grid::unit(c(0,1,0,1), "mm"))
  
  return(list(umap = umap/clonal_mixture + plot_layout(heights = c(1.45,1,0.1)),
              legend = umap_legend))
  
}

# function to run standard Seurat pre-processing and dim reduction with Immature cells
runHSCs <- function(seurat, ctypes = c("HSCs & MPPs"), ct_simpl = F, resolution = 0.8) {
  
  if(ct_simpl){
    
    W <- subset(seurat, ct %in% ctypes)
    
  }else{W <- subset(seurat, celltype %in% ctypes)}
  W <- NormalizeData(W)
  W <- ScaleData(W)
  W <- FindVariableFeatures(W)
  W <- RunPCA(W)
  W <- RunUMAP(W, reduction = "pca", dims = 1:10)
  W <- FindNeighbors(W, reduction = "pca", dims = 1:10)
  W <- FindClusters(W, resolution = resolution)
  #W$simple <- Idents(W)
  W
}

# function to plot expression of CD11c, CD49f
plot_markers <- function(seurat, leuk_prob = 0.9, healthy_prob = 0.9,
                         ct = "myeloid", cd11_thr = 0.3, cd49_thr = 1){
  
  # normalise RNA and ADT counts
  seurat <-NormalizeData(seurat)
  seurat <- NormalizeData(seurat, assay = "ADT", normalization.method = "CLR")
  
  # make data frame with RNA and antibody counts together with clonal assignments
  markers_data <- data.frame(cell_barcode = colnames(seurat),
                             CD34 = GetAssayData(seurat, slot = "data", 
                                                 assay = "ADT")["CD34.1",],
                             CD14 = GetAssayData(seurat, slot = "data", 
                                                 assay = "ADT")["CD14.1",],
                             CD49f = GetAssayData(seurat, slot = "data", 
                                                  assay = "ADT")["CD49f",],
                             CD11c = GetAssayData(seurat, slot = "data", 
                                                  assay = "RNA")["ITGAX",]) %>% 
    left_join(seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>% 
                dplyr::select(cell_barcode, ct, clone, leukemia_prob)) %>% 
    mutate(cancer = ifelse(clone == "Healthy", F, T)) %>% 
    filter(leukemia_prob > leuk_prob | leukemia_prob < (1-healthy_prob))
  
  if(ct == "myeloid"){
    
    markers_data <- markers_data %>% 
                      filter(!ct %in% c("T cells", "B cells", "NK cells", "Plasma cells", "Other")) 
      
  }else if(ct == "Immature"){
    
    markers_data <- markers_data %>% 
                      filter(ct == "Immature")
    
  }else if(ct == "CD34+"){
    
    markers_data <- markers_data %>% 
                      filter(CD34 > 1)
    
  }else if(ct == "CD14+"){
    
    markers_data <- markers_data %>% 
      filter(CD14 > 1)
    
  }
  
  # discretisize cell based on marker values
  markers_data <- markers_data %>% 
    mutate(gate = case_when(CD49f > cd49_thr & CD11c < cd11_thr ~ "CD49f+CD11c-",
                            CD49f < cd49_thr & CD11c > cd11_thr ~ "CD49f-CD11c+",
                            CD49f > cd49_thr & CD11c > cd11_thr ~ "CD49F+CD11c+",
                            CD49f < cd49_thr & CD11c < cd11_thr ~ "CD49F-CD11c-"))
  
  markers_fisher <- markers_data %>% filter(gate %in% c("CD49f+CD11c-", "CD49f-CD11c+")) %>% 
                      group_by(gate, cancer) %>% 
                      tally() %>% 
                      ungroup() %>% 
                      pivot_wider(id_cols = cancer, 
                                  names_from = gate, values_from = n) 
  
  # check whether there are cells in all gates
  if(which(is.na(markers_fisher)) %>% length() > 0){
    
    print("Not enough cells in all gates")
    
    return(markers_fisher)
    
  }else{markers_fisher <- fisher.test(markers_fisher)}
    
  # add p_value to dataframe
  markers_data <- markers_data %>% mutate(pval = paste0(pat, " (p_value = ",
                                                        markers_fisher$p.value, ")"))
  
  gates_plot <- ggplot(markers_data, aes(x = CD49f, y = CD11c, color = cancer)) +
    geom_point(size = 0.95, position = position_jitter(width=0.1,height=0.1)) +
    theme_classic() +
    scale_color_brewer(palette = "Set1") +
    geom_vline(xintercept = cd49_thr) +
    ggtitle(ct) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = cd11_thr) +
    facet_wrap(~pval)
  
  return(gates_plot)
  
}

# function to make boxplot for the expression of certain genes among groups
box_expression <- function(seurat, genes, barcodes = NULL, column = "status"){
  
  # if no barcodes are provide all cells in Seurat object are used
  if(is.null(barcodes)){
    
    bc <- colnames(seurat)
    
  }
  
  # normalise counts from Seurat object
  seurat <- NormalizeData(seurat)
  
  # make data frame with expression of genes
  data <- seurat@assays$RNA@data[genes,] %>% t() %>% as.data.frame() 
  
  colnames(data) <- genes
  
  # add clonal identities
  data <- data %>% mutate(cell_barcode = colnames(seurat)) %>% 
    left_join(seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>% 
                dplyr::select(cell_barcode,one_of(column))) %>% 
    filter(cell_barcode %in% bc) %>% 
    pivot_longer(cols = one_of(genes), names_to = "gene", values_to = "expression")
  
  # make plot
  boxplot <- ggplot(data, aes_string(x = column, y = "expression")) +
    geom_boxplot()+
    theme_classic()+
    ylab("Norm. expression") +
    facet_wrap(vars(gene))
  
  return(boxplot)
  
}

# I wrapped all the commands in a function to make it easier if several comparison want to be done
mast_dea <- function(barcodes, seurat, covariates, variable,
                     lrt, output, fdr_thr = 2, fold_thr = log2(1.5)){
  
  
  # get raw counts from the seurat object for the cells of interest
  counts <- seurat@assays$RNA@counts[,barcodes]
  
  # put gene names in a data frame
  fdata <- data.frame(row.names=rownames(counts), key = rownames(counts))
  
  # put data in a format compatible with MAST
  sca <- FromMatrix(as.matrix(log2(counts+1)), 
                    seurat@meta.data[colnames(counts),], fdata)
  
  # put the variable and covariates in the formula format for zlm 
  model <- paste0("~", paste(c(variable, covariates), collapse = "+"))
  
  # run the model with log10 of number of features and celltype as covariates
  zlmCond <- zlm(as.formula(model), sca, parallel = T)
  
  
  #then one can fiddle together the result, the output that MAST gives is not very user friendly...
  de_summ <- summary(zlmCond, doLRT= lrt)  
  
  # get p_values
  fdr_table <- subset(de_summ$datatable, component == "H", select= c("primerid","contrast","Pr(>Chisq)"))
  
  # adjust for multiple comparisons using FDR correction
  fdr_table$fdr <- p.adjust(fdr_table$`Pr(>Chisq)`, method = "fdr")
  
  # get fold change 
  fold_change <- subset(de_summ$datatable, component == "logFC", select = c("primerid","contrast","coef"))
  
  # merge adjusted p_values and fold changes
  final_table <- merge(fdr_table, fold_change, by = c("primerid","contrast")) %>% as.data.frame() %>% 
    mutate(log_fdr = -log10(fdr),
           log2_fold = coef) %>% 
    arrange(desc(log_fdr)) %>% 
    mutate(significant = ifelse((log_fdr > fdr_thr & abs(log2_fold) > fold_thr), TRUE, FALSE))
  
  # save final table as csv
  write_csv(final_table, file = output)
  
  return(list(model_obj = de_summ, sig_table = final_table))
  
}
