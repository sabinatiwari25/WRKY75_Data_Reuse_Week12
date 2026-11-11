source('scripts/setup_library.R')
dir.create('02_Validation', showWarnings = FALSE)
packages <- c('DESeq2','ggplot2','pheatmap','jsonlite','org.At.tair.db','GO.db',
              'AnnotationDbi','BiocManager','readxl','clusterProfiler','rmarkdown','knitr','renv')
installed <- vapply(packages, requireNamespace, logical(1), quietly=TRUE)
versions <- vapply(packages, function(p) if (installed[[p]])
  as.character(getNamespaceVersion(p)) else NA_character_, character(1))
write.csv(data.frame(package=packages, loadable=installed, version=versions),
          '02_Validation/software.csv', row.names=FALSE)
writeLines(sub('[[:blank:]]+$','',capture.output(sessionInfo())), '02_Validation/R_session_info.txt')
print(data.frame(package=packages, loadable=installed, version=versions))
