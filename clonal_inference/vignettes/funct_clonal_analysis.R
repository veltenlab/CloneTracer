# scripts with functions to carry out clonal analysis of patients
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(reticulate)
pd <- import("pandas")
library(circlize)
library(qpdf)
library(DiagrammeRsvg)
library(RColorBrewer)
library(ArchR)
library(rsvg)
library(DiagrammeR)

warriors_yellow <- "#FDBB30"
healthy_blue <- "#3b4ba7"
leuk_red <- "#ed1f26ff"
ind_color <- "#58595bff"
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

# function to add simple celltypes to Seurat object
add_simple_ct <- function(seurat){
  
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
  
  
  
  seurat$ct <- ifelse(seurat$ct %in% names(ct2simple), ct2simple[seurat$ct], "Other")
  
  return(seurat)
  
  
}

# function to make a dataframe with all information about a particular sample (after running the probabilistic model)
tidy_pickle <- function(pickle, clone_names, tree = "0"){
  
  data <- data.frame(cell_barcode = pickle$cell_barcode,
                     clone_n = apply(pickle$clonal_prob[[tree]], 1, which.max),
                     clonal_probability = apply(pickle$clonal_prob[[tree]], 1, max)) %>%  
    mutate(clone = factor(clone_names[clone_n], levels = clone_names),
           leukemia_prob = 1-pickle$clonal_prob[[tree]][,1], 
           status = ifelse(clone_n != 1, "cancer", "healthy")) %>% 
    arrange(clone_n)
  
  return(data)
  
}

# function to fix mito names
mito_rename <- function(string){
  
  pos <- gsub("X([0-9]+).+$", "\\1", string)
  
  ref <- gsub(".+\\.([A-Z]).+", "\\1", string)
  
  alt <- gsub(".+\\.([A-Z])$", "\\1", string)
  
  name <- paste0("mt:", pos, ref, ">", alt)
  
}

# function to get attributes for selected tree
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
  start_ind <- end_nodes <- rep("", nrow(tree_mat)-1)
  
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
      
      # get parent node
      node_parent <- tail(node_order[node_order %in% as.character(parents[[nxt_node]]+1)], n = 1)
      
      # get new mutations in the node
      node_name <- paste(colnames(tree_mat)[which(tree_mat[nxt_node,]-tree_mat[as.integer(node_parent),] == 1)], 
                         collapse = "\n")
      
      # add name of node and the index
      node_names[node] <- node_name
      node_order[node] <- nxt_node
      
      # get parent node
      start_nodes[node-1]  <- node_names[which(node_order == node_parent)] 
      
      end_nodes[node-1] <- node_name
      
    }
  }
  
  tree_mat <- tree_mat[as.integer(node_order),]
  rownames(tree_mat) <- node_names
  
  return(list("start_nodes" = start_nodes, "end_nodes" = end_nodes,
              "tree_mat" = tree_mat, "clone_names" = node_names[order(node_order)]))
  
}

  

# change the tree plotting function to include new discrete colours
# function to plot all clonal hierarchies selected by the model
plot_trees <- function(pickle, tree = "all", clone_cols = F, fontype = "bold", 
                       cols = c("#3b4ba7", rev(ArchRPalettes$calm)),
                       fontcols = c(rep("#FFFFFF", 3), rep("#000000", 2), "#FFFFFF", "#000000", "#FFFFFF", rep("#000000", 2))){
  
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
        fixed_size = FALSE,
        fontname = "Helvetica",
        fontsize = 11, 
        fontcolor = '#FFFFFF',
        fontype = fontype,
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
    
    # if(clone_cols){
    #   
    #   combs <- combn(nrow(tree_att$tree_mat), 1, function(x)replace(numeric(nrow(tree_att$tree_mat)),x,1))
    #   colors <- rev(apply(combs, 1, getColors))
    #   
    # }else{colors <- "black"}
    colors <- cols[1:length(rownames(tree_att$tree_mat))]
    
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
      fontcolor = fontcols[1:length(rownames(tree_att$tree_mat))],
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
  
  # create output directory if not present
  suppressWarnings(dir.create(outdir))
  
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

# I have to change the default function to order cells by leukemia probability
make_heatmap <- function(pickle, seurat, pat, tree, cols = c("#3b4ba7", rev(ArchRPalettes$calm)),
                         cnv_type = NULL, cnv_pos = NULL, cnv_matrix = NULL, middle_point = 0.5,
                         cvg = F, cvg_max = 10, cvg_mid = 5, clust_rows = F,
                         prob_type = "clones", celltype = F, ct_column = NULL){
  
  # get order of nodes in the tree and clone names
  tree_att <- get_tree_attr(pickle, tree)
  
  # get clonal probs and name them
  probs <- pickle$clonal_prob[[tree]]
  colnames(probs) <- tree_att$clone_names
  probs <- probs[,rownames(tree_att$tree_mat)]
  rownames(probs) <- pickle$cell_barcode
  pickle$clonal_prob[[tree]] <- probs
  clone_names <- tree_att$clone_names
  probs <- probs[, rev(1:length(clone_names))] 
  
  # get discrete clone colours
  clone_cols_discr <- rev(cols[1:length(clone_names)])
  
  # clone names in the order of clonal probability matrix
  clones_order <- rownames(tree_att$tree_mat)
  
  # put data into a dataframe for plotting 
  table <- tidy_pickle(pickle,
                        clone_names = clones_order,
                        tree = tree) %>% 
                        mutate(patient = pat) 
  
  # celltype info ---------------------------------------------------------------------------------------------------
  # add celltype identity if specified
  if(celltype){
    
    # make sure same cell barcodes are present
    bc <- intersect(colnames(seurat), table$cell_barcode)
    
    # add celltype column from Seurat object
    table <- table %>% filter(cell_barcode %in% bc) %>% 
              left_join(seurat@meta.data %>% rownames_to_column(var = "cell_barcode") %>% 
                          dplyr::select(cell_barcode, one_of(ct_column)))
    
    # create pallete for celltypes (up to 29 groups)
    ct_cols <- c("#0C727C", "#7E1416", ArchRPalettes$summerNight, rev(ArchRPalettes$stallion))[1:length(unique(table[[ct_column]]))]
    
    if(prob_type == "clones"){
      
      # order cells based on clonal probabilities
      table <- table %>% 
                group_by(clone) %>% 
                arrange(desc(clonal_probability), .by_group = T)
      
      # make barplot of clonal probabilities
      barplot <- anno_barplot(probs[table$cell_barcode,],
                              gp = gpar(fill = clone_cols_discr, col = NA), 
                              bar_width = 1, height = unit(1.9, "cm"), which = "column", 
                              border = T, axis_param = list(at = c(0,0.25,0.5,0.75,1)))
      
      # create object with the heatmap column annotations
      col_ann <- HeatmapAnnotation("Clonal\nprobability" = barplot,
                                   "Celltype" = table[[ct_column]],
                                   show_legend = c(F, T),
                                   col = list("Celltype" = setNames(ct_cols, unique(table[[ct_column]]))),
                                   simple_anno_size = unit(0.7, "cm"),
                                   annotation_legend_param = list(Celltype = list(nrow = 2, title_position = "topcenter")))
      
      # create matrix with VAF to plot in the heatmap
      heatmap_matrix <- pickle$M/(pickle$N+pickle$M)
      
      # rename mitochondrial mutation names
      mutnames <- sapply(pickle$mutations_matrix, function(i){
        
        ifelse(grepl("\\.", i), mito_rename(i), i)
        
      })
      
      # Unit scale trisomy reads and rename cnvs
      if (length(cnv_type) > 0){
        
        for(i in 1:length(cnv_type)){
          
          max <- max(pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]])
          
          scaled_ratio <- pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]]/max
          
          heatmap_matrix[,cnv_pos[i]] <- scaled_ratio
          
        }
      }
      
      names(mutnames) <- NULL
      
      colnames(heatmap_matrix) <- mutnames
      rownames(heatmap_matrix) <- pickle$cell_barcode
      
      heatmap_matrix <- heatmap_matrix[table$cell_barcode,] %>% t()
      
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
                         height = unit(1.4*nrow(heatmap_matrix), "cm"),
                         show_column_names = F,
                         use_raster = F,
                         na_col = dropout_gray,
                         column_title = "Cells", 
                         column_title_side = "bottom") 
      
      # add legend with clone names and also dropout status
      lgd_list <- list(Legend(labels = "Dropout", legend_gp = gpar(fill = dropout_gray), title = "Status"),
                       Legend(labels = rownames(tree_att$tree_mat),
                              legend_gp = gpar(fill = rev(clone_cols_discr)),
                              title = "Clone"))

    }else if(prob_type == "cancer"){
      
      # order cells based on cancer probability
      table <- table %>% 
                arrange(leukemia_prob)
      
      # make heatmap with cancer probability
      probs_mat <- table$leukemia_prob %>% as.matrix() %>% t()
      
      # heatmap object
      heat_probs <- Heatmap(probs_mat, 
                            heatmap_legend_param = list(title = "Cancer\nProbability", legend_direction = "vertical",
                                                        legend_width = unit(3, "cm"),
                                                        title_position = "topcenter"),
                            col = colorRamp2(c(0,0.5,1), c(healthy_blue, ind_color, leuk_red)),
                            show_heatmap_legend = T,
                            cluster_columns = F,
                            cluster_rows = F,
                            show_column_dend = F,
                            show_row_dend = F, 
                            show_row_names = F,
                            height = unit(1.9, "cm"),
                            show_column_names = F,
                            use_raster = F,
                            na_col = dropout_gray)
      
      # make heatmap with VAF 
      # create object with the heatmap column annotations
      col_ann <- HeatmapAnnotation("Celltype" = table[[ct_column]],
                                   show_legend = c(T),
                                   col = list("Celltype" = setNames(ct_cols, unique(table[[ct_column]]))),
                                   simple_anno_size = unit(0.7, "cm"),
                                   annotation_legend_param = list(Celltype = list(nrow = 2, title_position = "topcenter")))
      
      # create matrix with VAF to plot in the heatmap
      heatmap_matrix <- pickle$M/(pickle$N+pickle$M)
      
      # rename mitochondrial mutation names
      mutnames <- sapply(pickle$mutations_matrix, function(i){
        
        ifelse(grepl("\\.", i), mito_rename(i), i)
        
      })
      
      # Unit scale trisomy reads and rename cnvs
      if (length(cnv_type) > 0){
        
        for(i in 1:length(cnv_type)){
          
          max <- max(pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]])
          
          scaled_ratio <- pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]]/max
          
          heatmap_matrix[,cnv_pos[i]] <- scaled_ratio
          
        }
      }
      
      names(mutnames) <- NULL
      
      colnames(heatmap_matrix) <- mutnames
      rownames(heatmap_matrix) <- pickle$cell_barcode
      
      heatmap_matrix <- heatmap_matrix[table$cell_barcode,] %>% t()
      
      # make heatmap in which the order of the rows is given by the clonal information
      heat_vaf <- Heatmap(heatmap_matrix,
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
                         height = unit(1.4*nrow(heatmap_matrix), "cm"),
                         show_column_names = F,
                         use_raster = F,
                         na_col = dropout_gray,
                         column_title = "Cells", 
                         column_title_side = "bottom") 
      
      # add legend with clone names and also dropout status
      lgd_list <- list(Legend(labels = "Dropout", legend_gp = gpar(fill = dropout_gray), title = "Status"))
      
      heatmap = heat_probs %v% heat_vaf
      
    }else{
      
      # order cells based on clonal probabilities
      table <- table %>% 
        group_by(clone) %>% 
        arrange(desc(clonal_probability), .by_group = T)
      
      # make barplot of clonal probabilities
      barplot <- anno_barplot(probs[table$cell_barcode,],
                              gp = gpar(fill = clone_cols_discr, col = NA), 
                              bar_width = 1, height = unit(1.9, "cm"), which = "column", 
                              border = T, axis_param = list(at = c(0,0.25,0.5,0.75,1)))
      
      # make heatmap with cancer probability
      probs_mat <- table$leukemia_prob %>% as.matrix() %>% t()
      
      # heatmap object
      heat_probs <- Heatmap(probs_mat, 
                            heatmap_legend_param = list(title = "Cancer\nProbability", legend_direction = "vertical",
                                                        legend_width = unit(3, "cm"),
                                                        title_position = "topcenter"),
                            col = colorRamp2(c(0,0.5,1), c(healthy_blue, ind_color, leuk_red)),
                            show_heatmap_legend = T,
                            cluster_columns = F,
                            cluster_rows = F,
                            show_column_dend = F,
                            show_row_dend = F, 
                            show_row_names = F,
                            height = unit(1.9, "cm"),
                            show_column_names = F,
                            use_raster = F,
                            na_col = dropout_gray)
      
      # make heatmap with VAF 
      # create object with the heatmap column annotations
      col_ann <- HeatmapAnnotation("Clonal\nprobability" = barplot,
                                   "Celltype" = table[[ct_column]],
                                   show_legend = c(F, T),
                                   col = list("Celltype" = setNames(ct_cols, unique(table[[ct_column]]))),
                                   simple_anno_size = unit(0.7, "cm"),
                                   annotation_legend_param = list(Celltype = list(nrow = 2, title_position = "topcenter")))
      
      # create matrix with VAF to plot in the heatmap
      heatmap_matrix <- pickle$M/(pickle$N+pickle$M)
      
      # rename mitochondrial mutation names
      mutnames <- sapply(pickle$mutations_matrix, function(i){
        
        ifelse(grepl("\\.", i), mito_rename(i), i)
        
      })
      
      # Unit scale trisomy reads and rename cnvs
      if (length(cnv_type) > 0){
        
        for(i in 1:length(cnv_type)){
          
          max <- max(pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]])
          
          scaled_ratio <- pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]]/max
          
          heatmap_matrix[,cnv_pos[i]] <- scaled_ratio
          
        }
      }
      
      names(mutnames) <- NULL
      
      colnames(heatmap_matrix) <- mutnames
      rownames(heatmap_matrix) <- pickle$cell_barcode
      
      heatmap_matrix <- heatmap_matrix[table$cell_barcode,] %>% t()
      
      # make heatmap in which the order of the rows is given by the clonal information
      heat_vaf <- Heatmap(heatmap_matrix,
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
                         height = unit(1.4*nrow(heatmap_matrix), "cm"),
                         show_column_names = F,
                         use_raster = F,
                         na_col = dropout_gray,
                         column_title = "Cells", 
                         column_title_side = "bottom") 
      
      # add legend with clone names and also dropout status
      lgd_list <- list(Legend(labels = "Dropout", legend_gp = gpar(fill = dropout_gray), title = "Status"),
                       Legend(labels = rownames(tree_att$tree_mat),
                              legend_gp = gpar(fill = rev(clone_cols_discr)),
                              title = "Clone"))
      
      
      heatmap = heat_probs %v% heat_vaf
      
    }
    
  # no celltype info -------------------------------------------------------------------------------------------------------
  }else{
    
    if(prob_type == "clones"){
      
      # order cells based on clonal probabilities
      table <- table %>% 
        group_by(clone) %>% 
        arrange(desc(clonal_probability), .by_group = T)
      
      # make barplot of clonal probabilities
      barplot <- anno_barplot(probs[table$cell_barcode,],
                              gp = gpar(fill = clone_cols_discr, col = NA), 
                              bar_width = 1, height = unit(1.9, "cm"), which = "column", 
                              border = T, axis_param = list(at = c(0,0.25,0.5,0.75,1)))
      
      # create object with the heatmap column annotations
      col_ann <- HeatmapAnnotation("Clonal\nprobability" = barplot,
                                   show_legend = c(F),
                                   simple_anno_size = unit(0.7, "cm"))
      
      # create matrix with VAF to plot in the heatmap
      heatmap_matrix <- pickle$M/(pickle$N+pickle$M)
      
      # rename mitochondrial mutation names
      mutnames <- sapply(pickle$mutations_matrix, function(i){
        
        ifelse(grepl("\\.", i), mito_rename(i), i)
        
      })
      
      # Unit scale trisomy reads and rename cnvs
      if (length(cnv_type) > 0){
        
        for(i in 1:length(cnv_type)){
          
          max <- max(pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]])
          
          scaled_ratio <- pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]]/max
          
          heatmap_matrix[,cnv_pos[i]] <- scaled_ratio
          
        }
      }
      
      names(mutnames) <- NULL
      
      colnames(heatmap_matrix) <- mutnames
      rownames(heatmap_matrix) <- pickle$cell_barcode
      
      heatmap_matrix <- heatmap_matrix[table$cell_barcode,] %>% t()
      
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
                         height = unit(1.4*nrow(heatmap_matrix), "cm"),
                         show_column_names = F,
                         use_raster = F,
                         na_col = dropout_gray,
                         column_title = "Cells", 
                         column_title_side = "bottom") 
      
      # add legend with clone names and also dropout status
      lgd_list <- list(Legend(labels = "Dropout", legend_gp = gpar(fill = dropout_gray), title = "Status"),
                       Legend(labels = rownames(tree_att$tree_mat),
                              legend_gp = gpar(fill = rev(clone_cols_discr)),
                              title = "Clone"))
      
    }else if(prob_type == "cancer"){
      
      # order cells based on cancer probability
      table <- table %>% 
        arrange(leukemia_prob)
      
      # make heatmap with cancer probability
      probs_mat <- table$leukemia_prob %>% as.matrix() %>% t()
      
      # heatmap object
      heat_probs <- Heatmap(probs_mat, 
                            heatmap_legend_param = list(title = "Cancer\nProbability", legend_direction = "vertical",
                                                        legend_width = unit(3, "cm"),
                                                        title_position = "topcenter"),
                            col = colorRamp2(c(0,0.5,1), c(healthy_blue, ind_color, leuk_red)),
                            show_heatmap_legend = T,
                            cluster_columns = F,
                            cluster_rows = F,
                            show_column_dend = F,
                            show_row_dend = F, 
                            show_row_names = F,
                            height = unit(1.9, "cm"),
                            show_column_names = F,
                            use_raster = F,
                            na_col = dropout_gray)
      
      # make heatmap with VAFs
      # create matrix with VAF to plot in the heatmap
      heatmap_matrix <- pickle$M/(pickle$N+pickle$M)
      
      # rename mitochondrial mutation names
      mutnames <- sapply(pickle$mutations_matrix, function(i){
        
        ifelse(grepl("\\.", i), mito_rename(i), i)
        
      })
      
      # Unit scale trisomy reads and rename cnvs
      if (length(cnv_type) > 0){
        
        for(i in 1:length(cnv_type)){
          
          max <- max(pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]])
          
          scaled_ratio <- pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]]/max
          
          heatmap_matrix[,cnv_pos[i]] <- scaled_ratio
          
        }
      }
      
      names(mutnames) <- NULL
      
      colnames(heatmap_matrix) <- mutnames
      rownames(heatmap_matrix) <- pickle$cell_barcode
      
      heatmap_matrix <- heatmap_matrix[table$cell_barcode,] %>% t()
      
      # make heatmap in which the order of the rows is given by the clonal information
      heat_vaf <- Heatmap(heatmap_matrix,
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
                         height = unit(1.4*nrow(heatmap_matrix), "cm"),
                         show_column_names = F,
                         use_raster = F,
                         na_col = dropout_gray,
                         column_title = "Cells", 
                         column_title_side = "bottom") 
      
      # add legend with clone names and also dropout status
      lgd_list <- list(Legend(labels = "Dropout", legend_gp = gpar(fill = dropout_gray), title = "Status"))
      
      heatmap = heat_probs %v% heat_vaf
      
    }else{
      
      # order cells based on clonal probabilities
      table <- table %>% 
        group_by(clone) %>% 
        arrange(desc(clonal_probability), .by_group = T)
      
      # make barplot of clonal probabilities
      barplot <- anno_barplot(probs[table$cell_barcode,],
                              gp = gpar(fill = clone_cols_discr, col = NA), 
                              bar_width = 1, height = unit(1.9, "cm"), which = "column", 
                              border = T, axis_param = list(at = c(0,0.25,0.5,0.75,1)))
      
      # make heatmap with cancer probability
      probs_mat <- table$leukemia_prob %>% as.matrix() %>% t()
      
      # heatmap object
      heat_probs <- Heatmap(probs_mat, 
                            heatmap_legend_param = list(title = "Cancer\nProbability", legend_direction = "vertical",
                                                        legend_width = unit(3, "cm"),
                                                        title_position = "topcenter"),
                            col = colorRamp2(c(0,0.5,1), c(healthy_blue, ind_color, leuk_red)),
                            show_heatmap_legend = T,
                            cluster_columns = F,
                            cluster_rows = F,
                            show_column_dend = F,
                            show_row_dend = F, 
                            show_row_names = F,
                            height = unit(1.9, "cm"),
                            show_column_names = F,
                            use_raster = F,
                            na_col = dropout_gray)
      
      # make heatmap with VAF 
      # create object with the heatmap column annotations
      col_ann <- HeatmapAnnotation("Clonal\nprobability" = barplot,
                                   show_legend = c(F),
                                   simple_anno_size = unit(0.7, "cm"))
      
      # create matrix with VAF to plot in the heatmap
      heatmap_matrix <- pickle$M/(pickle$N+pickle$M)
      
      # rename mitochondrial mutation names
      mutnames <- sapply(pickle$mutations_matrix, function(i){
        
        ifelse(grepl("\\.", i), mito_rename(i), i)
        
      })
      
      # Unit scale trisomy reads and rename cnvs
      if (length(cnv_type) > 0){
        
        for(i in 1:length(cnv_type)){
          
          max <- max(pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]])
          
          scaled_ratio <- pickle$M[,cnv_pos[i]]/pickle$N[,cnv_pos[i]]/max
          
          heatmap_matrix[,cnv_pos[i]] <- scaled_ratio
          
        }
      }
      
      names(mutnames) <- NULL
      
      colnames(heatmap_matrix) <- mutnames
      rownames(heatmap_matrix) <- pickle$cell_barcode
      
      heatmap_matrix <- heatmap_matrix[table$cell_barcode,] %>% t()
      
      # make heatmap in which the order of the rows is given by the clonal information
      heat_vaf <- Heatmap(heatmap_matrix,
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
                         height = unit(1.4*nrow(heatmap_matrix), "cm"),
                         show_column_names = F,
                         use_raster = F,
                         na_col = dropout_gray,
                         column_title = "Cells", 
                         column_title_side = "bottom") 
      
      # add legend with clone names and also dropout status
      lgd_list <- list(Legend(labels = "Dropout", legend_gp = gpar(fill = dropout_gray), title = "Status"),
                       Legend(labels = rownames(tree_att$tree_mat),
                              legend_gp = gpar(fill = rev(clone_cols_discr)),
                              title = "Clone"))
      
      
      heatmap = heat_probs %v% heat_vaf
    }
    
  }
  
  # coverage 
  # plot heatmap showing total coverage underneath
  if(cvg){
    
    # create matrix with VAF to plot in the heatmap
    heatmap_matrix <- pickle$N+pickle$M
    
    mutnames <- sapply(pickle$mutations_matrix, function(i){
      
      ifelse(grepl("\\.", i), mito_rename(i), i)
      
    })
    
    # Rename cnvs (Otherwise it throws an error)
    if (length(cnv_type) > 0){
      
      for(i in 1:length(cnv_type)){
        
        mutnames[cnv_pos[i]] <- paste0(cnv_type[i], "\n", mutnames[cnv_pos[i]])
        
      }
    }
    
    names(mutnames) <- NULL
    
    colnames(heatmap_matrix) <- mutnames
    rownames(heatmap_matrix) <- pickle$cell_barcode
    
    heatmap_matrix <- heatmap_matrix[table$cell_barcode,] %>% t()
    
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
                           height = unit(1*nrow(heatmap_matrix), "cm"),
                           show_column_names = F,
                           use_raster = F,
                           na_col = dropout_gray,
                           column_title = "Cells", 
                           column_title_side = "bottom") 
    
    ht_list = heatmap %v% heatmap_cvg
    
    #draw(ht_list, annotation_legend_list = lgd_list, merge_legend = TRUE)
    
    return(list(plot = ht_list, legend = lgd_list))
    
  }else{
    
    return(list(plot = heatmap, legend = lgd_list))  
    
  }
}

# function to make UMAP with clones
umap_clones <- function(pickle, seurat, tree, name, post_thr = 0.8, reduction = "umap"){
  
  # get order of nodes in the tree and clone names
  tree_att <- get_tree_attr(pickle, tree)
  
  clone_names <- tree_att$clone_names
  
  # get clonal probs and name them
  probs <- pickle$clonal_prob[[tree]]
  colnames(probs) <- tree_att$clone_names
  probs <- probs[,rownames(tree_att$tree_mat)]
  rownames(probs) <- pickle$cell_barcode
  pickle$clonal_prob[[tree]] <- probs
  
  # clone names in the order of clonal probability matrix
  clones_order <- rownames(tree_att$tree_mat)
  
  # get coordinates
  coords <- as.data.frame(Embeddings(seurat, reduction = reduction)) %>% 
    rownames_to_column(var = "cell_barcode")
  
  colnames(coords) <- c("cell_barcode", paste(reduction, 1, sep = "_"), paste(reduction, 2, sep = "_"))
  
  # a cell is assigned to a clone if the posterior probability is higher than the threshold (default 0.8)
  data <- tidy_pickle(pickle,
                      clone_names = clones_order,
                      tree = tree) %>%  
                    mutate(clone = ifelse(clonal_probability > post_thr, as.character(clone), "Indeterminate")) %>% 
                    left_join(coords) %>% 
                    filter(cell_barcode %in% intersect(colnames(seurat), pickle$cell_barcode))
  
  # set colours to clones following the inferred clonal hierarchy
  cols <- setNames(c(c("#3b4ba7", rev(ArchRPalettes$calm))[1:length(clones_order)], ind_color),
                   c(clones_order, "Indeterminate"))
  
  # UMAP with reference as background
  umap <- ggplot(data,
                 aes_string(x = paste(reduction, 1, sep = "_"), y = paste(reduction, 2, sep = "_"),
                            color = "clone")) +
                  geom_point(size = 0.65) +
                  theme_classic() +
                  scale_color_manual(values = cols) +
                  #coord_cartesian(xlim = c(-13, 15), ylim = c(-12, 12.5)) +
                  theme(axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.line = element_blank(),
                        panel.border = element_rect(color = "black", fill = NA),
                        axis.ticks = element_blank(),
                        legend.position = "bottom",
                        legend.title = element_blank(),
                        plot.margin=grid::unit(c(0.75,0,0,0), "mm"),
                        strip.text = element_text(size = 13)) +
                  guides(color = guide_legend(override.aes = list(size = 4)))
  
    
    return(umap)
 
}

# function to make UMAP with leukemia probability
umap_leuk <- function(pickle, seurat, tree, reduction = "umap"){
  
    # get order of nodes in the tree and clone names
    tree_att <- get_tree_attr(pickle, tree)
    
    clone_names <- tree_att$clone_names
    
    # get clonal probs and name them
    probs <- pickle$clonal_prob[[tree]]
    colnames(probs) <- tree_att$clone_names
    probs <- probs[,rownames(tree_att$tree_mat)]
    rownames(probs) <- pickle$cell_barcode
    pickle$clonal_prob[[tree]] <- probs
    
    # clone names in the order of clonal probability matrix
    clones_order <- rownames(tree_att$tree_mat)
  
    # get coordinates
    coords <- as.data.frame(Embeddings(seurat, reduction = reduction)) %>% 
                rownames_to_column(var = "cell_barcode")
    
    colnames(coords) <- c("cell_barcode", paste(reduction, 1, sep = "_"), paste(reduction, 2, sep = "_"))
    
    # put data in data.frame and add coordinates for dim. reduction plot
    data <- tidy_pickle(pickle,
                        clone_names = clones_order,
                        tree = tree) %>% 
            left_join(coords) %>% 
            filter(cell_barcode %in% intersect(colnames(seurat), pickle$cell_barcode))
              
    
    # UMAP with reference as background
    umap <- ggplot(data,
                   aes_string(x = paste(reduction, 1, sep = "_"), y = paste(reduction, 2, sep = "_"), 
                              color = "leukemia_prob")) +
                  geom_point(size = 0.65) +
                  theme_classic() +
                  scale_color_gradient2(low = healthy_blue, mid = ind_color, high = leuk_red, midpoint = 0.5) +
                  #coord_cartesian(xlim = c(-13, 15), ylim = c(-12, 12.5)) +
                  theme(axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.line = element_blank(),
                        panel.border = element_rect(color = "black", fill = NA),
                        axis.ticks = element_blank(),
                        legend.position = "bottom",
                        plot.margin=grid::unit(c(0.75,0,0,0), "mm"),
                        strip.text = element_text(size = 13)) +
                  guides(color = guide_colorbar(title = "Cancer probability", title.position = "top", 
                                                barwidth = unit(4, "cm")))
    

      
      return(umap)

}

# function to add clonal metadata to Seurat
add_meta <- function(seurat, pickle, tree){
  
  tree_att <- get_tree_attr(pickle, tree)
  
  clone_names <- gsub("\n", "-", tree_att$clone_names)
  
  clone_to_names <- setNames(clone_names, 1:length(clone_names))
  
  probs <- pickle$clonal_prob[[tree]]
  
  bc <- intersect(colnames(seurat), pickle$cell_barcode)
  
  meta <- seurat@meta.data %>% 
    rownames_to_column(var = "cell_barcode") %>% 
    left_join(data.frame(cell_barcode = pickle$cell_barcode,
                         clone = apply(probs, 1, which.max),
                         leukemia_prob = 1-probs[,1],
                         clonal_probability = apply(probs, 1, max))) %>%
    filter(cell_barcode %in% bc) %>% 
    mutate(clone = clone_to_names[clone],
           status = ifelse(clone == "Healthy", "healthy", "cancer")) %>% 
    column_to_rownames(var = "cell_barcode") %>% 
    dplyr::rename(cancer_probability = leukemia_prob)
  
  seurat@meta.data <- meta
  
  return(seurat)
  
}


