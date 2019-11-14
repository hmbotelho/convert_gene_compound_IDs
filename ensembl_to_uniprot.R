# ensembl_to_uniprot_chr
#
# Converts human ensembl gene symbols to Uniprot IDs. Parallel processing available
# 
# Hugo Botelho & André Falcão
# v0.1
# 14 November 2019
#
# Input: 
#	* ens: character vector, with human ensembl gene IDs
#	* showProgress: logical, print progress to console?
#	* parallelize: logical, use multi-processors?
# Output: character vector, with Uniprot IDs
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


ensembl_to_uniprot_chr <- function(ens, showProgress = FALSE, parallelize = FALSE){
        
	# Make sure a mart is available
    if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
	
    ens <- as.character(ens)
    
    if(detectCores() <= 2 | !parallelize)
    {
        # Single-core processing
        output <- lapply(ens, function(ENS){
            
            if(showProgress) print(paste0("Converting '", ENS, "' to Uniprot ID."), quote = FALSE)
                
            if(exists("gene_table")) rm("gene_table")
            tryCatch({
                gene_table <- getBM(attributes = c('ensembl_gene_id', 'uniprotswissprot'),
                                    filters = 'ensembl_gene_id', values = ENS, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(ensembl_gene_id = character(0),
                                             uniprotswissprot = character(0))
                }
            })
            
            if(nrow(gene_table) == 0){
                result <- ""
            } else{
                result <- gene_table[1,"uniprotswissprot"]
            }
            
            result
            
        })
        
    } else
    {
        # Multi-core processing
        cluster <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cluster, varlist=c("ens", "myMart", "getBM"), envir=environment())
        
        output <- parLapply(cluster, ens, function(ENS){
            
            if(exists("gene_table")) rm("gene_table")
            tryCatch({
                gene_table <- getBM(attributes = c('ensembl_gene_id', 'uniprotswissprot'),
                                    filters = 'ensembl_gene_id', values = ENS, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(ensembl_gene_id = character(0),
                                             uniprotswissprot = character(0))
                }
            })
            
            if(nrow(gene_table) == 0){
                result <- ""
            } else{
                result <- gene_table[1,"uniprotswissprot"]
            }
            
            result
            
        })
        stopCluster(cluster)
    }
    
    # Return character vector
    output <- do.call("rbind", output)
    as.character(output)
}