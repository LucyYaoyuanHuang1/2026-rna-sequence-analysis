#setup----
  #set working directory, for clarity, following Qi
  setwd("C:/Users/Lyhuang/Downloads/Data")
  #load packages, magrittr for piping, dplyr and tibble for data manip, chm for heatmap, circ for heatmap colors
  library("tidyverse")
  library("ComplexHeatmap")
  library("circlize")
  library("ClassComparison")
  library("effectsize")
  library("magrittr")

#read in and prepare data----
  sample_data <- 
    read.csv("Samples.csv") %>% 
    mutate(Sample = str_replace_all(Sample,"-",".")) %>% 
    as_tibble()
  #read in pool data, set the first column as rownames, as matrix, apply log(x + 1)
  pool_data <- 
    read.csv("pool.tpm_Gene.csv") %>%  
    column_to_rownames(var = "X") %>% 
    mutate(across(where(is.double), ~add(.x, 1))) %>% 
    mutate(across(where(is.double), log2)) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Sample")
  a_data <- pool_data %>% 
    semi_join(sample_data %>% 
    filter(TimePoint == "48H"), by = join_by(Sample)) %>% 
    column_to_rownames(var = "Sample") %>% 
    t()
  a_classes <- 
    sample_data %>% 
    filter(TimePoint == "48H") %>% 
    select(Group) %>% pull(Group) %>% 
    factor()
  b_data <- 
    pool_data %>% 
    semi_join(sample_data %>% 
    filter(TimePoint == "72H"), by = join_by(Sample)) %>% 
    column_to_rownames(var = "Sample") %>% 
    t()
  b_classes <- 
    sample_data %>% 
    filter(TimePoint == "72H") %>% 
    select(Group) %>% pull(Group)%>% 
    factor()
  c_data <- 
    pool_data %>% 
    semi_join(sample_data %>% 
    filter(Treatment == "TTFields"), by = join_by(Sample)) %>% 
    column_to_rownames(var = "Sample") %>% 
    t()
  c_classes <- 
    sample_data %>% 
    filter(Treatment == "TTFields") %>% 
    select(Group) %>% pull(Group)%>% 
    factor()
  Pvalcut<-c(0.00001,0.0001,0.001,0.005,0.01,0.05,seq(0.1,1,by=0.05))

#multi t-testing 48H----
  a_mttest <- MultiTtest(a_data, a_classes, na.rm = FALSE)
  a_bum <- Bum(a_mttest@p.values)
  hist(a_bum, main = "48H")
  a_significant_genes <- a_data[selectSignificant(a_bum, alpha = .05, by = "FDR"),] %>% na.omit()
  
  #multi t-testing 72H----
  b_mttest <- MultiTtest(b_data, b_classes, na.rm = TRUE)
  b_bum <- Bum(b_mttest@p.values)
  hist(b_bum, main = "72H")
  b_significant_genes <- b_data[selectSignificant(b_bum, alpha = .05,by = "FDR"),] %>% na.omit()
  
  #multi t-testing TTFields----
  c_mttest <- MultiTtest(c_data, c_classes, na.rm = FALSE)
  c_bum <- Bum(c_mttest@p.values)
  hist(c_bum, main = "TTFields")
  c_significant_genes <- c_data[selectSignificant(c_bum, alpha = .05,by = "FDR"),] %>% na.omit()


#heatmaps----
  
  #general and reference----
  
  #for ease of constant plot creation, use head as test data, and standardize for legibility of graph (+ and - values for diff colors)
  #small_pool_data = a_significant_genes[1:10,] %>% standardize()
  
  #set colors to be between red and blue and between -2 and 2
  color_function <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  #make a heatmap, colors using previous color_function, row dendrogram width
  #made much smaller, row names made much smaller, cluster_columns between groups (preserving overall order) by 
  #sample_data$group(possible sources of error: using sample_data for grouping, read documentation)

    column_heatmap_annotation_48H = HeatmapAnnotation(
      group = a_classes
    )
    heatmap_48h_control_ttfields <- Heatmap(
      a_significant_genes %>% standardize(),
      top_annotation = column_heatmap_annotation_48H,
      name = "Group", 
      col = color_function, 
      row_dend_width = unit(.1, "cm"), 
      row_names_gp = gpar(fontsize = 2),
    )
    heatmap_48h_control_ttfields
    
  
#heatmap notes----
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
  #general ComplexHeatmap info

