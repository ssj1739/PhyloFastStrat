# Installing all dependencies

install.packages('devtools')

devtools::install_github("arendsee/phylostratr")

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
  "tidyverse"
  )

BiocManager::install(package_list)
