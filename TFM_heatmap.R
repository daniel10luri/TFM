install.packages("pheatmap")
library(pheatmap)
library(dplyr)
library(tidyr)


# Cargar archivos
argannot_raw <- read_tsv("argannot_final.txt")
ncbi_raw <- read_tsv("ncbi_final.txt")
resfinder_raw <- read_tsv("resfinder_final.txt")
card_raw <- read_tsv("card_final.txt")

# Función corregida para preparar los datos
prepare_data <- function(df) {
  df %>%
    select(`#FILE`, `GENE`) %>% 
    distinct() %>%         
    pivot_wider(names_from = GENE, values_from = GENE, values_fill = list(GENE = "")) %>%  
    mutate(across(-`#FILE`, ~ ifelse(. != "", 1, 0))) %>%  
    mutate(across(-`#FILE`, as.integer)) %>%  
    as.data.frame()  
}
# Preparar los datos para la base de datos 
argannot_heatmap <- prepare_data(argannot)
card_heatmap <- prepare_data(card)
ncbi_heatmap <- prepare_data(ncbi)
resfinder_heatmap <- prepare_data(resfinder)

# Unir los 4 dataframes y seleccionar solo FILE y GENE
alldb <- bind_rows(argannot, ncbi, resfinder, card) %>%
  select(`#FILE`, GENE) %>%
  mutate(GENE = tolower(GENE)) %>%   
  distinct(`#FILE`, GENE, .keep_all = TRUE) 
View(alldb)


prepare_all_data <- function(df) {
  df %>%
  
    pivot_wider(names_from = GENE, values_from = GENE, values_fill = list(GENE = "")) %>%  
    mutate(across(-FILE, ~ ifelse(. != "", 1, 0))) %>%  
    mutate(across(-FILE, as.integer)) %>% 
    as.data.frame() 
}

# Aplicar la función al dataframe
alldb_heatmap <- prepare_data(alldb)

# Convertimos la primera columna en nombres de fila para evitar el problema
rownames(alldb_heatmap) <- alldb_heatmap$`#FILE`
alldb_heatmap <- alldb_heatmap[, -1]  

# Calcular la prevalencia de cada gen (suma de las columnas)
gene_prevalence_alldb <- colSums(alldb_heatmap)
ordered_genes_alldb <- names(sort(gene_prevalence_alldb, decreasing = TRUE))
alldb_heatmap <- alldb_heatmap[, ordered_genes_alldb]

#Crear el gráfico
pheatmap(as.matrix(alldb_heatmap), 
         clustering_distance_rows = "manhattan",  
         clustering_distance_cols = "manhattan" ,clustering_method = "complete",
         color = colorRampPalette(c("beige", "lightsalmon"))(50),
         fontsize_row = 8,   
         fontsize_col = 6, treeheight_row = 100,  
         treeheight_col = 100, 
        legend = TRUE, border_color = NA)     


















