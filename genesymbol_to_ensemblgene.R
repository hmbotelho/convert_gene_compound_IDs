# genesymbol_to_ensemblgene_df
#
# Converts human gene symbols to Ensembl gene IDs
# 
# Hugo Botelho & André Falcão
# v0.1
# 14 November 2019
#
# Input: 
#	* genes: character vector, with human gene symbols
# Output: data frame, with the input and corresponding Ensembl gene IDs
#
# Dependencies: biomaRt


if(!("biomaRt" %in% installed.packages())){
	source("https://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
}
if(!("biomaRt" %in% loadedNamespaces())) library(biomaRt)
if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")    # This can take a while

if(!("parallel" %in% installed.packages())) install.packages("parallel")
library(parallel)

genesymbol_to_ensemblgene_df <- function(genes) {
    
    # Make sure a mart is available
    if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
    
    results_table<-getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'), 
                    filters = 'hgnc_symbol', values = genes, mart = myMart)
    results_table<-results_table[nchar(results_table$ensembl_gene_id)>0,]
    results_table
}







# genesymbol_to_ensemblgene_chr
#
# Converts human gene symbols to Ensembl gene IDs
# 
# Hugo Botelho & André Falcão
# v0.1
# 14 November 2019
#
# Input: 
#	* genes: character vector, with human gene symbols
#   * showProgress: logical, print progress to console?
#	* parallelize: logical, use multi-processors?
# Output: character vector, with the input and corresponding Ensembl gene IDs
#
# Dependencies: biomaRt, parallel

genesymbol_to_ensemblgene_chr <- function(genes, showProgress = FALSE, parallelize = FALSE){
        
	# Make sure a mart is available
    if(!exists("myMart")) myMart <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
	
    genes <- as.character(genes)
    
    if(detectCores() <= 2 | !parallelize)
    {
        # Single-core processing
        output <- lapply(genes, function(GENE){
            
            if(showProgress) print(paste0("Converting '", GENE, "' to Ensembl ID."), quote = FALSE)
                
            if(exists("gene_table")) rm("gene_table")
            tryCatch({
                gene_table <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                                    filters = 'hgnc_symbol', values = GENE, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(hgnc_symbol = character(0),
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
        clusterExport(cluster, varlist=c("genes", "myMart", "getBM"), envir=environment())
        
        output <- parLapply(cluster, genes, function(GENE){
            
            if(exists("gene_table")) rm("gene_table")
            tryCatch({
                gene_table <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                                    filters = 'hgnc_symbol', values = GENE, mart = myMart)
                gene_table <- gene_table[gene_table[,2]!="",]
            },
            error = function(e){},
            finally = {
                if(!exists("gene_table")){
                    gene_table <- data.frame(hgnc_symbol = character(0),
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