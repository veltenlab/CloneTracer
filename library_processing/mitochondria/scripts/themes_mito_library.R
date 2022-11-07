# Here themes and other ggplot specifications are specified for bulk ATAC plots as well as colours

library(RColorBrewer)
library(tidyverse)

google_red <- "#DB4437"
google_blue <- "#4285F4"  
google_yellow <- "#F4B400"
google_green <- "#0F9D58"
mitochondria_purple <- brewer.pal(5, "Set1")[4]
other_colour <- brewer.pal(5, "Set1")[5]
orange_set1 <- brewer.pal(8, "Dark2")[2]
golden <- brewer.pal(8, "Dark2")[7]
salmon <- brewer.pal(3, "Accent")[3]
dropout_gray <- "#d3d3d3"


col_cvg_comparison <- c(google_blue, google_red, google_green, google_yellow, "black", "gray")


# theme for alignment plots in rmarkdown
theme_plots_alignment <-  list(theme_classic() +
                                     theme(axis.line.x = element_blank(),
                                           plot.title = element_text(size = 20, hjust = 0.5, vjust = 0.5),
                                           legend.title = element_blank(),
                                           legend.text = element_text(size = 11),
                                           legend.key.size = unit(20, "pt"),
                                           axis.ticks.x = element_blank(),
                                           axis.title.y = element_text(size = 14),
                                           axis.title.x = element_blank(), 
                                           axis.text = element_text(size = 11, 
                                                                    colour = "black", hjust = 0.5)))


# theme to save alignment plots
theme_plots_alignment_pdf <-  list(theme_classic() +
                                 theme(axis.line.x = element_blank(),
                                       plot.title = element_text(size = 24, hjust = 0.5, vjust = 0.5),
                                       legend.title = element_blank(),
                                       legend.text = element_text(size = 13),
                                       legend.key.size = unit(26, "pt"),
                                       axis.ticks.x = element_blank(),
                                       axis.title.y = element_text(size = 16),
                                       axis.title.x = element_blank(), 
                                       axis.text = element_text(size = 14, 
                                                                colour = "black", hjust = 0.5)))


# theme for raw coverage plots in rmarkdown
theme_raw_coverage <-  list(theme_classic(),
                            theme(axis.line = element_blank(),
                                  legend.position = "bottom",
                                  legend.title = element_blank(),
                                  legend.text = element_text(size = 11),
                                  legend.key.size = unit(20, "pt"),
                                  axis.ticks.x = element_blank(),
                                  axis.title = element_text(size = 14),
                                  axis.text.y = element_text(size = 11, 
                                                             colour = "black", hjust = 0.5)),
                            scale_x_continuous(expand = c(0,50), 
                                               breaks = c(1, seq(1000, 15000, 1000))),
                            scale_y_continuous(trans = "log10",
                                               breaks = c(1,10,100, 1000),
                                               labels = c("1x", "10x", "100x", "1000x")),
                            geom_hline(yintercept = c(1, 10, 100, 1000), 
                                       linetype = "dashed",
                                       alpha = 0.3)) 

# theme to save raw coverage plots 
theme_raw_coverage_pdf <-  list(theme_classic(),
                                theme(axis.line = element_blank(),
                                      legend.position = "bottom",
                                      legend.title = element_blank(),
                                      legend.text = element_text(size = 13),
                                      legend.key.size = unit(26, "pt"),
                                      axis.ticks.x = element_blank(),
                                      axis.title = element_text(size = 16),
                                      axis.text.y = element_text(size = 14, 
                                                                 colour = "black", hjust = 0.5)),
                                scale_x_continuous(expand = c(0,50), 
                                                   breaks = c(1, seq(1000, 15000, 1000))),
                                scale_y_continuous(trans = "log10",
                                                   breaks = c(1,10,100, 1000),
                                                   labels = c("1x", "10x", "100x", "1000x")),
                                geom_hline(yintercept = c(1, 10, 100, 1000), 
                                           linetype = "dashed",
                                           alpha = 0.3)) 



# theme for raw coverage plots in rmarkdown
theme_norm_coverage <-  list(theme_classic(),
                              theme(axis.line = element_blank(),
                                    legend.position = "bottom",
                                    legend.title = element_blank(),
                                    legend.text = element_text(size = 11),
                                    legend.key.size = unit(20, "pt"),
                                    axis.ticks.x = element_blank(),
                                    axis.title = element_text(size = 14),
                                    axis.text.y = element_text(size = 11, 
                                                               colour = "black", hjust = 0.5)),
                              scale_x_continuous(expand = c(0,50), 
                                                breaks = c(1, seq(1000, 15000, 1000))),
                              scale_y_continuous(trans = "log10"),
                              geom_hline(yintercept = c(1e-2, 1e-4, 1e-6), 
                                        linetype = "dashed",
                                        alpha = 0.3)) 


# theme to save normalised coverage plots 
theme_norm_coverage_pdf <-  list(theme_classic(),
                                 theme(axis.line = element_blank(),
                                       legend.position = "bottom",
                                       legend.title = element_blank(),
                                       legend.text = element_text(size = 13),
                                       legend.key.size = unit(26, "pt"),
                                       axis.ticks.x = element_blank(),
                                       axis.title = element_text(size = 16),
                                       axis.text.y = element_text(size = 14, 
                                                                  colour = "black", hjust = 0.5)),
                                 scale_x_continuous(expand = c(0,50), 
                                                    breaks = c(1, seq(1000, 15000, 1000))),
                                 scale_y_continuous(trans = "log10"),
                                 geom_hline(yintercept = c(1e-2, 1e-4, 1e-6), 
                                            linetype = "dashed",
                                            alpha = 0.3)) 


# theme for raw coverage plots in rmarkdown
theme_umi_coverage <-  list(theme_classic(),
                            theme(axis.line = element_blank(),
                                  legend.position = "bottom",
                                  legend.title = element_blank(),
                                  legend.text = element_text(size = 11),
                                  legend.key.size = unit(20, "pt"),
                                  axis.ticks.x = element_blank(),
                                  axis.title = element_text(size = 14),
                                  axis.text.y = element_text(size = 11, 
                                                             colour = "black", hjust = 0.5)),
                            scale_x_continuous(expand = c(0,50), 
                                               breaks = c(1, seq(1000, 15000, 1000))),
                            scale_y_continuous(trans = "log10",
                                               breaks = c(1,3,10,50),
                                               labels = c("1x", "3x","10x","50x")),
                            geom_hline(yintercept = c(1, 3, 10, 50), 
                                       linetype = "dashed",
                                       alpha = 0.3)) 

# theme to save raw coverage plots 
theme_umi_coverage_pdf <-  list(theme_classic(),
                                theme(axis.line = element_blank(),
                                      legend.position = "bottom",
                                      legend.title = element_blank(),
                                      legend.text = element_text(size = 13),
                                      legend.key.size = unit(26, "pt"),
                                      axis.ticks.x = element_blank(),
                                      axis.title = element_text(size = 16),
                                      axis.text.y = element_text(size = 14, 
                                                                 colour = "black", hjust = 0.5)),
                                scale_x_continuous(expand = c(0,50), 
                                                   breaks = c(1, seq(1000, 15000, 1000))),
                                scale_y_continuous(trans = "log10",
                                                   breaks = c(1,3,10,50),
                                                   labels = c("1x", "3x","10x","50x")),
                                geom_hline(yintercept = c(1, 3, 10, 50), 
                                           linetype = "dashed",
                                           alpha = 0.3)) 

# theme for base quality plot in rmarkdown
theme_base_qual <-  list(theme_classic(),
                            theme(axis.line = element_blank(),
                                  legend.position = "bottom",
                                  legend.title = element_blank(),
                                  legend.text = element_text(size = 11),
                                  legend.key.size = unit(20, "pt"),
                                  axis.ticks.x = element_blank(),
                                  axis.title = element_text(size = 14),
                                  axis.text.y = element_text(size = 11, 
                                                             colour = "black", hjust = 0.5)),
                            scale_x_continuous(expand = c(0,50), 
                                               breaks = c(1, seq(1000, 15000, 1000))),
                            scale_y_continuous(breaks = c(0,10, 20, 30, 40)),
                            geom_hline(yintercept = c(1, 10, 20, 30, 40), 
                                       linetype = "dashed",
                                       alpha = 0.3)) 

# theme for raw coverage plots in rmarkdown
theme_base_qual_pdf <-  list(theme_classic(),
                             theme(axis.line = element_blank(),
                                   legend.position = "bottom",
                                   legend.title = element_blank(),
                                   legend.text = element_text(size = 13),
                                   legend.key.size = unit(26, "pt"),
                                   axis.ticks.x = element_blank(),
                                   axis.title = element_text(size = 16),
                                   axis.text.y = element_text(size = 14, 
                                                              colour = "black", hjust = 0.5)),
                         scale_x_continuous(expand = c(0,50), 
                                            breaks = c(1, seq(1000, 15000, 1000))),
                         scale_y_continuous(breaks = c(0,10, 20, 30, 40)),
                         geom_hline(yintercept = c(1, 10, 20, 30, 40), 
                                    linetype = "dashed",
                                    alpha = 0.3)) 
