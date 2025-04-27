# Installing all dependencies

install.packages('devtools')

devtools::install_github("arendsee/phylostratr", upgrade = T)

install.packages("BiocManager")

package_list <- c(
  "shiny",
  "shinycssloaders",
  "queryup",
  "readxl",
  "org.Hs.eg.db",
  "devtools",
  "tidyverse",
  "topGO",
  "GO.db",
  "fgsea",
  "AnnotationDbi",
  "pbapply",
  "STRINGdb",
  "quarto",
  "DT",
  "tidyverse",
  "googlesheets4",
  "ggtext",
  "shinythemes"
  )

BiocManager::install(package_list, update = T, ask = F)
devtools::install_github("alserglab/fgsea", upgrade = F)

# uniprotkb_organism_id_9606_2025_03_15_tsv <- read_delim("~/Documents/JainAnalytics/Kolabtree/Proposal1/uniprotkb_organism_id_9606_2025_03_15.tsv.gz", 
#                                                         delim = "\t", escape_double = FALSE, 
#                                                         trim_ws = TRUE)

