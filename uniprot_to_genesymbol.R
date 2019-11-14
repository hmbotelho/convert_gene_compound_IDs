# uniprot_to_genesymbol_chr
#
# Converts Uniprot IDs to human gene symbols. Parallel processing available
# 
# Hugo Botelho & André Falcão
# v0.1
# 14 November 2019
#
# Input: 
#	* up: character vector, with Uniprot IDs
#	* showProgress: logical, print progress to console?
#	* parallelize: logical, use multi-processors?
# Output: character vector, with human gene symbols
#
# Dependencies: biomaRt, parallel


if(!("biomaRt" %in% installed.packages())){
	source("https://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
}
if(!("biomaRt" %in% loadedNamespaces())) library(biomaRt)
if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")    # This can take a while

if(!("parallel" %in% installed.packages())) install.packages("parallel")
library(parallel)




uniprot_to_genesymbol_chr <- function(up, showProgress = FALSE, parallelize = FALSE){
    	
	# Make sure a mart is available
    if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
    
    up <- as.character(up)
    
    if(detectCores() <= 2 | !parallelize)
    {
        # Single-core processing
        output <- lapply(up, function(UP){
            
            if(showProgress) print(paste0("Converting '", UP, "' to gene symbol."), quote = FALSE)
                
            if(exists("gene_table")) rm("gene_table")
            tryCatch({
                gene_table <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol'),
                                    filters = 'uniprotswissprot', values = UP, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(uniprotswissprot = character(0),
                                             hgnc_symbol = character(0))
                }
            })
            
            if(nrow(gene_table) == 0){
                result <- ""
            } else{
                result <- gene_table[1,"hgnc_symbol"]
            }
            
            result
            
        })
        
    } else
    {
        # Multi-core processing
        cluster <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cluster, varlist=c("up", "myMart", "getBM"), envir=environment())
        
        output <- parLapply(cluster, up, function(UP){
            
            if(exists("gene_table")) rm("gene_table")
            tryCatch({
                gene_table <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol'),
                                    filters = 'uniprotswissprot', values = UP, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(uniprotswissprot = character(0),
                                             hgnc_symbol = character(0))
                }
            })
            
            if(nrow(gene_table) == 0){
                result <- ""
            } else{
                result <- gene_table[1,"hgnc_symbol"]
            }
            
            result
            
        })
        stopCluster(cluster)
    }
    
    # Return character vector
    output <- do.call("rbind", output)
    as.character(output)
}