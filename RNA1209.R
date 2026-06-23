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

#read in sample data----
sample_data <-  read.csv("Samples.csv")
#read in pool data, set the first column as rownames, as matrix, apply log(x + 1) (addition first to avoid -Inf, might be cause of no negatives?)
pool_data = read.csv("pool.tpm_Gene.csv") %>%  column_to_rownames(var = "X") %>% as.matrix() %>% add(1) %>% log2()


#multi t-testing48H----
  A_classes <- sample_data %>% as_tibble() %>% filter(Group == "Control 48H" | Group == "TTFields 48H") %>% select(Group) %>% pull(Group)%>% factor()
  A_mttest <- MultiTtest(pool_data, A_classes)
  summary(A_mttest)
  A_bum <- Bum(A_mttest@p.values)
  Pvalcut<-c(0.00001,0.0001,0.001,0.005,0.01,0.05,seq(0.1,1,by=0.05))
  summary(A_bum, Pvalcut)
  hist(A_bum, main = "48H")
  A_signifigant_genes <- selectSignificant(A_bum, alpha = .05,by = "FDR")
  
  #multi t-testing72H----
  B_classes <- sample_data %>% as_tibble() %>% filter(Group == "Control 48H" | Group == "TTFields 48H") %>% select(Group) %>% pull(Group)%>% factor()
  B_mttest <- MultiTtest(pool_data, B_classes)
  summary(B_mttest)
  B_bum <- Bum(B_mttest@p.values)
  Pvalcut<-c(0.00001,0.0001,0.001,0.005,0.01,0.05,seq(0.1,1,by=0.05))
  summary(B_bum, Pvalcut)
  hist(B_bum, main = "72H")
  B_signifigant_genes <- selectSignificant(B_bum, alpha = .05,by = "FDR")


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
    column_heatmap_annotation_48h = HeatmapAnnotation(
      group = sample_data$Group, 
      col = list(group = c(
      "TTFields 48H" = "red",
      "Control 48H" = "blue"),
      height = unit(200, "cm")))
  
  
    heatmap_48h_control_ttfields <- Heatmap(
      small_pool_data,
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

