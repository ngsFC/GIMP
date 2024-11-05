.onLoad <- function(libname, pkgname) {
    # List of Bioconductor packages package depends on
    bioc_packages <- c("IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylationEPICanno.ilm10b4.hg19","IlluminaHumanMethylationEPICv2anno.20a1.hg38") 
    # Check if BiocManager is available; install it if not
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }

    # Install missing Bioconductor dependencies
    for (pkg in bioc_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            BiocManager::install(pkg, ask = FALSE)
        }
    }
    
    packageStartupMessage("Thank you for using the GIMP package. Dont forget to cite us!")
}
