# few lines to help make sure we do have everything we need
`%not in%` <- Negate(`%in%`)

needed_packages <-
  c("remotes", "devtools", "BiocManager", "knitr", "rmarkdown",
    "tidyverse", "dplyr", "tidyr",
    "dada2", "vegan", "picante",
    "ggplot2", "patchwork", "scales", "rstatix",
    "exactRankTests", "compositions", "nlme",
    "gtools", "Rcpp", "igraph")

# a helper abbreviation
`%not in%` <- Negate(`%in%`)

for (pkg in needed_packages) {
  if (pkg %not in% rownames(installed.packages())) {
    print(paste("Trying to install", pkg))
    install.packages(pkg)
    if ((pkg %not in% rownames(installed.packages()))) {
      msg <- paste("ERROR: Unsuccessful!", pkg, "not installed!",
                   "Check the log and try installing the package manually.")
      stop(msg)
    } 
  }
  library(pkg, character.only = TRUE)
  ifelse(pkg %in% loadedNamespaces(), 
         print(paste("Successful!", pkg, "loaded.")),
         print(paste("ERROR: Unsuccessful!", pkg, 
                     "not loaded. Check error msg.")))
}

# additional packages from BiocManager
needed_packages <- c("ALDEx2","DECIPHER","DESeq2","edgeR","phyloseq")
for (pkg in needed_packages) {
  if (pkg %not in% rownames(installed.packages())) {
    print(paste("Trying to install", pkg))
    BiocManager::install(pkg)
    if ((pkg %not in% rownames(installed.packages()))) {
      msg <- paste("ERROR: Unsuccessful!", pkg, "not installed!",
                   "Check the log and try installing the package manually.")
      stop(msg)
    } 
  }
  library(pkg, character.only = TRUE)
  ifelse(pkg %in% loadedNamespaces(), 
         print(paste("Successful!", pkg, "loaded.")),
         print(paste("ERROR: Unsuccessful!", pkg, 
                     "not loaded. Check error msg.")))
}

# additional packages from github
needed_packages_remotes <- 
  c("KlausVigo/phangorn", "kassambara/ggpubr", "thomasp85/patchwork",
    "YinchengChen23/PreLectR", "zdk123/SpiecEasi")



for (pkg in needed_packages_remotes) {
  pkg_name <- basename(pkg)
  if (pkg_name %not in% rownames(installed.packages())) {
    print(paste("Trying to install", pkg_name))
    remotes::install_github(pkg)
    if (pkg_name %not in% rownames(installed.packages())) {
      msg <- paste("ERROR: Unsuccessful!", pkg, "not installed!",
                   "Check the log and try installing the package manually.")
      stop(msg)
    } 
  }
  library(pkg_name, character.only = TRUE)
  ifelse(pkg_name %in% loadedNamespaces(), 
         print(paste("Successful!", pkg, "loaded.")),
         print(paste("ERROR: Unsuccessful!", pkg, 
                     "not loaded. Check error msg.")))
}
