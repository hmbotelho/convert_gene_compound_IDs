# uniprot_to_ensembl_chr
#
# Converts Uniprot IDs to Ensembl Gene IDs. Parallel processing available
# 
# Hugo Botelho & André Falcão
# v0.1
# 14 November 2019
#
# Input: 
#	* up: character vector, with Uniprot IDs
#	* showProgress: logical, print progress to console?
#	* parallelize: logical, use multi-processors?
# Output: character vector, with Ensembl gene IDs
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



uniprot_to_ensembl_chr <- function(up, showProgress = FALSE, parallelize = FALSE){
	
	# Make sure a mart is available
    if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
    
    up <- as.character(up)
    
    if(detectCores() <= 2 | !parallelize)
    {
        # Single-core processing
        output <- lapply(up, function(UP){
            
            if(showProgress) print(paste0("Converting '", UP, "' to Ensembl ID."), quote = FALSE)
                
            if(exists("gene_table")) rm("gene_table")
            tryCatch({
                gene_table <- getBM(attributes = c('uniprotswissprot', 'ensembl_gene_id'),
                                    filters = 'uniprotswissprot', values = UP, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(uniprotswissprot = character(0),
                                             ensembl_gene_id = character(0))
                }
            })
            
            if(nrow(gene_table) == 0){
                result <- ""
            } else{
                result <- gene_table[1,"ensembl_gene_id"]
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
                gene_table <- getBM(attributes = c('uniprotswissprot', 'ensembl_gene_id'),
                                    filters = 'uniprotswissprot', values = UP, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(uniprotswissprot = character(0),
                                             ensembl_gene_id = character(0))
                }
            })
            
            if(nrow(gene_table) == 0){
                result <- ""
            } else{
                result <- gene_table[1,"ensembl_gene_id"]
            }
            
            result
            
        })
        stopCluster(cluster)
    }
    
    # Return character vector
    output <- do.call("rbind", output)
    as.character(output)
}