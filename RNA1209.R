#setup----
  #set working directory, for clarity, following Qi
  setwd("C:/Users/Lyhuang/Downloads/Data")
  library("tidyverse")
  library("ComplexHeatmap")
  library("circlize")
  library("ClassComparison")
  library("effectsize")
  library("magrittr")
  library("qqplotr")
  library("patchwork")

#read in and prepare data----
  sample_data <- 
    read.csv("Samples.csv") %>% 
    mutate(Sample = str_replace_all(Sample,"-",".")) %>% 
    as_tibble()
  #read in pool data, set the first column as rownames, as matrix, apply log(x + 1)
  pool_data <- 
    read.csv("pool.tpm_Gene.csv") %>%  
    column_to_rownames(var = "X") %>%
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Sample")
  
  log_pool_data <- 
    read.csv("pool.tpm_Gene.csv") %>%  
    column_to_rownames(var = "X") %>% 
    mutate(across(where(is.double), ~add(.x, 1))) %>% 
    mutate(across(where(is.double), log2)) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Sample")
  
#taking a look at base data distribution----
  long_pool = pool_data %>% pivot_longer(cols = !Sample, names_to = "gene", values_to = "values") %>% filter(values >0) %>% slice(1:400)
  log_long_pool = log_pool_data %>% pivot_longer(cols = !Sample, names_to = "gene", values_to = "values") %>% filter(values >0) %>% slice(1:400)
  no_zero_log = pool_data %>% filter(rowSums(across(everything()) != 0) > 0)
  
  hist <- ggplot(long_pool, mapping = aes(x = values)) + geom_histogram(bins = 100, fill = "green", alpha = .5) + 
    xlab("ggplot(long_pool, mapping = aes(x = values)) + geom_histogram(bins = 100, fill = 'green', alpha = .5)") + 
    labs(title = "Histogram, 1st 400 non-log adjusted values")
  qqplot <- ggplot(long_pool, mapping = aes(sample = values)) + geom_qq_band(fill = 'green', alpha = .5) + stat_qq_point() + stat_qq_line() +
      xlab("ggplot(long_pool, mapping = aes(sample = values)) + geom_qq_band(fill = 'green', alpha = .5) + stat_qq_point() + stat_qq_line()") + 
      labs(title = "QQ Plot, 1st 400 non-log adjusted values")
  log_hist <- ggplot(log_long_pool, mapping = aes(x = values)) + geom_histogram(bins = 100, fill = 'blue', alpha = .5) + 
    xlab("ggplot(log_long_pool, mapping = aes(x = values)) + geom_histogram(bins = 100, fill = 'blue', alpha = .5)") +
    labs(title = "Histogram, 1st 400 log adjusted values")
  log_qqplot <- ggplot(log_long_pool, mapping = aes(sample = values)) + geom_qq_band(fill = 'blue', alpha = .5) + stat_qq_point() + stat_qq_line() + 
    xlab("ggplot(log_long_pool, mapping = aes(sample = values)) + geom_qq_band(fill = 'blue', alpha = .5) + stat_qq_point() + stat_qq_line()") + 
    labs(title = "QQ Plot, 1st 400 log adjusted values")

  normality_figure <- (hist | qqplot) / (log_hist | log_qqplot)
  
  hist
  qqplot
  log_hist
  log_qqplot
  normality_figure
  
  shapiro.test(long_pool$values)
  shapiro.test(log_long_pool$values)
  
#organizing data for ease of analysis----
  make_classes <- function(grouping, attribute){
    sample_data %>% 
      filter({{grouping}} == attribute) %>% 
      select(Group) %>% 
      pull(Group) %>% 
      factor() %>% 
      return()
  }
  a_classes <- make_classes(TimePoint, "48H")
  b_classes <- make_classes(TimePoint, "72H")
  c_classes <- make_classes(Treatment, "TTFields")
  
  make_data <- function(grouping, attribute) {
    pool_data %>% 
      semi_join(sample_data %>% 
        filter({{grouping}} == attribute), by = join_by(Sample)) %>% 
      column_to_rownames(var = "Sample") %>% 
      t()
  }
  a_data <- make_data(TimePoint, "48H")
  b_data <- make_data(TimePoint, "72H")
  c_data <- make_data(Treatment, "TTFields")
#multi t-testing 48H----
  mttest <- function(data, classes) {
    temp_mttest <- MultiTtest(data, classes, na.rm = FALSE)
    a_bum <- Bum(a_mttest@p.values)
    hist(a_bum, main = "48H")
    a_cutoff = cutoffSignificant(a_bum, alpha = .05, by = "FDR")
    a_significant_genes <- a_data[selectSignificant(a_bum, by = "FDR"),] %>% na.omit()
  }
  
  a_mttest <- MultiTtest(a_data, a_classes, na.rm = FALSE)
  a_bum <- Bum(a_mttest@p.values)
  hist(a_bum, main = "48H")
  a_cutoff = cutoffSignificant(a_bum, alpha = .05, by = "FDR")
  a_significant_genes <- a_data[selectSignificant(a_bum, by = "FDR"),] %>% na.omit()
  
  #multi t-testing 72H
  b_mttest <- MultiTtest(b_data, b_classes, na.rm = TRUE)
  b_bum <- Bum(b_mttest@p.values)
  hist(b_bum, main = "72H")
  b_significant_genes <- b_data[selectSignificant(b_bum, alpha = .01, by = "FDR"),] %>% na.omit()
  
  #multi t-testing TTFields
  c_mttest <- MultiTtest(c_data, c_classes, na.rm = FALSE)
  c_bum <- Bum(c_mttest@p.values)
  hist(c_bum, main = "TTFields")
  c_significant_genes <- c_data[selectSignificant(c_bum, alpha = .01,by = "FDR"),] %>% na.omit()

#Fold changes----
  geo_mean <- function(){
    log_pool_data %>% 
      full_join(sample_data %>% 
                  select(Group, Sample), by = join_by(Sample)) %>%
      group_by(Group) %>% 
      summarize_if(is.numeric, mean) %>% 
      return()
  }
  geo_means <- geo_mean()
  
  fold_change <- function (geo_mean_data, group_1, group_2){
    geo_mean_data %>% 
      pivot_longer(cols = !Group, names_to = "Gene", values_to = "data") %>% 
      group_by(Gene) %>% 
      mutate(foch = ({{group_1}} - {{group_2}})) %>% 
      return()
  }
#todo subtract group1 from group2
  fc = fold_change(geo_means, `Control 48H`, `TTFields 48H`)
  
  
#Volcano maps----
#heatmaps----
  #set colors to be between red and blue and between -2 and 2
  color_function <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  #make a heatmap, colors using previous color_function, row dendrogram width
  #made much smaller, row names made much smaller, cluster_columns between groups (preserving overall order) by 
  #sample_data$group(possible sources of error: using sample_data for grouping, read documentation)

    column_heatmap_annotation_48H = HeatmapAnnotation(
      Group = a_classes,
      col = list(Group = c("Control 48H" = 3, "TTFields 48H" = 6))
    )
    heatmap_48h_control_ttfields <- Heatmap(
      a_significant_genes  %>% standardize(),
      top_annotation = column_heatmap_annotation_48H,
      name = "Group", 
      col = color_function, 
      row_dend_width = unit(.1, "cm"), 
      row_names_gp = gpar(fontsize = 2),
    )
    heatmap_48h_control_ttfields