# Path to data storage 
external_path <- "/Users/juhui/Desktop/research/data/"

# Package Installments 
packages <- c("dplyr", "tidyr", "haven", "dplyr", "ggplot2", "fixest")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}