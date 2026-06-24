#setup----
  #set working directory, for clarity, following Qi
  setwd("C:/Users/Lyhuang/Downloads/Data")
  #load packages, magrittr for piping, dplyr and tibble for data manip, chm for heatmap, circ for heatmap colors
  library("magrittr")
  library("dplyr")
  library("tibble")
  library("ComplexHeatmap")
  library("circlize")
  library("ClassComparison")
  library("effectsize")
  library("stringr")

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
#  summary(a_mttest)
  a_bum <- Bum(a_mttest@p.values)
  hist(a_bum, main = "48H")
  a_signifigant_genes <- a_data[selectSignificant(a_bum, alpha = .05,by = "FDR"),] %>% na.omit()
  
  #multi t-testing 72H----
  #ask if na.rm should be false
  b_mttest <- MultiTtest(b_data, b_classes, na.rm = FALSE)
  summary(b_mttest)
  b_bum <- Bum(b_mttest@p.values)
  hist(b_bum, main = "72H")
  b_signifigant_genes <- b_data[selectSignificant(b_bum, alpha = .05,by = "FDR"),] %>% na.omit()
  
  #multi t-testing TTFields----
  c_mttest <- MultiTtest(c_data, c_classes, na.rm = FALSE)
  summary(c_mttest)
  c_bum <- Bum(c_mttest@p.values)
  hist(c_bum, main = "TTFields")
  c_signifigant_genes <- c_data[selectSignificant(c_bum, alpha = .05,by = "FDR"),] %>% na.omit()


#heatmaps----
  
  #general and reference----
  
  #for ease of constant plot creation, use head as test data, and standardize for legibility of graph (+ and - values for diff colors)
  small_pool_data = pool_data[1:10,] %>% standardize()
  
  #set colors to be between red and blue and between -2 and 2
  color_function <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  #figure out what this line does? comes straight from example. perhaps turnes [-3, 3] into sequence of colors, why?
  color_function(seq(-3, 3))
  
  
  #make a heatmap with small_pool_data, legend name testing, colors using previous color_function, row dendrogram width
  #made much smaller, row names made much smaller, cluster_columns between groups (preserving overall order) by 
  #sample_data$group(possible sources of error: using sample_data for grouping, read documentation)
  column_heatmap_annotation = HeatmapAnnotation(
    group = sample_data$Group, 
    col = list(group = c(
    "TTFields 48H" = "red",
    "Control 48H" = "yellow",
    "TTFields 72H" = "green",
    "Control 72H" = "blue")),
    height = unit(200, "cm"))
  
  
  heatmap_control_ttfields <- Heatmap(small_pool_data, 
    name = "Group", 
    col = color_function, 
    row_dend_width = unit(.1, "cm"), 
    row_names_gp = gpar(fontsize = 3), 
    cluster_columns = cluster_within_group(small_pool_data, sample_data$Treatment),
    top_annotation = column_heatmap_annotation)
  heatmap_control_ttfields
  #48H----
    forty_eight_hours = sample_data %>% filter(TimePoint == "48H")
    column_heatmap_annotation_48h = HeatmapAnnotation(
      group = sample_data$Group, 
      col = list(group = c(
      "TTFields 48H" = "red",
      "Control 48H" = "blue"),
      height = unit(200, "cm")))
  
    heatmap_48h_control_ttfields <- Heatmap(
      small_pool_data %>% as_tibble() %>% semi_join(sample_data, by = ),
      name = "Group", 
      col = color_function, 
      row_dend_width = unit(.1, "cm"), 
      row_names_gp = gpar(fontsize = 3), 
      cluster_columns = cluster_within_group(small_pool_data, sample_data$Group),
      top_annotation = column_heatmap_annotation_48h)
    
    heatmap_48h_control_ttfields
  


  
#heatmap notes----
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
  #general ComplexHeatmap info

