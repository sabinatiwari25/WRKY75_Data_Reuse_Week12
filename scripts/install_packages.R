source('scripts/setup_library.R')
dir.create('.Rlib',showWarnings=FALSE)
.libPaths(c(normalizePath('.Rlib'),.libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=600)
if (!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',lib='.Rlib')
cran <- c('ggplot2','pheatmap','jsonlite')
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,lib='.Rlib',Ncpus=2)
bio <- c('DESeq2','AnnotationDbi','org.At.tair.db','GO.db')
missing <- bio[!vapply(bio,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) BiocManager::install(missing,version='3.19',lib='.Rlib',ask=FALSE,update=FALSE,Ncpus=2)
stopifnot(all(vapply(c(cran,bio),requireNamespace,logical(1),quietly=TRUE)))
dir.create('02_Validation',showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),'02_Validation/install_session.txt')
