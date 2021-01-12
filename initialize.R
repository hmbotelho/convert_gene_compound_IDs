# Convert gene and compound IDs
#
# Initialize dependencies
#
# Hugo Botelho
# v0.2
# 12 January 2020


source("https://raw.githubusercontent.com/hmbotelho/generic_scripts/main/R/setupPackages.R")

# biomaRt
if(!("biomaRt" %in% installed.packages())){
    
    r_version_full <- paste0(R.version$major, ".", R.version$minor)
    if(compareVersion(r_version_full, "3.5") >= 0){
        # Install in R >= 3.5
        install.packages("BiocManager")
        BiocManager::install("biomaRt")
    }else{
        # Install in R < 3.5
        source("https://bioconductor.org/biocLite.R")
        biocLite("biomaRt")
    }
}
setupPackages("biomaRt")
setupPackages("rentrez")


if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")    # This can take a while


