# geneid_to_genesymbol_chr
#
# Converts NCBI gene IDs to official gene symbols.
# 
# Hugo Botelho
# v0.1
# 14 November 2019
#
# Input: 
#	* geneids: vector (character or numeric), with NCBI gene IDs
#	* showProgress: logical, print progress to console?
# Output: character vector, with gene symbols
#
# Dependencies: rentrez


if(!("rentrez" %in% installed.packages())) install.packages("rentrez")
library(rentrez)


geneid_to_genesymbol_chr <- function(geneids, showProgress = FALSE){
    
    #geneids <- 1080
        
    results <- sapply(geneids, function(geneid){
        
        if(showProgress) print(paste0("Converting '", geneid, "' to Uniprot ID"), quote=FALSE)
        
        # Check that the refseq identifier exists and is unambiguous
        preliminary_info <- entrez_search(db="gene", term = geneid)
        if(preliminary_info$count != 1) return(NULL)
        
        
        # Get the most updated record for this gene
        summary1 <- entrez_summary(db="gene", id=geneid)
        if("replacedby" %in% names(summary1)){
            geneid_updated <- summary1$replacedby
        } else{
            geneid_updated <- geneid
        }
        
        
        # Get all Uniprot IDs from the record
        ncbi_data <- entrez_fetch(db="gene", id=geneid_updated, rettype="raw")
        
        regex_genesymbol <- ".*Official Symbol: (.*?) .*"
        genesymbol_loc0 <- gregexpr(regex_genesymbol, ncbi_data)
        genesymbol_loc  <- data.frame(start = genesymbol_loc0[[1]],
                                      end   = genesymbol_loc0[[1]] + attr(genesymbol_loc0[[1]], which="match.length") - 1)
        genesymbols <- apply(genesymbol_loc, 1, function(x){
            sub(regex_genesymbol, "\\1", substring(ncbi_data, x["start"], x["end"]))
        })
        
        
        # Select the most likely ID (i.e. the most abundant)
        id_table <- table(genesymbols)
        genesymbol <- names(which(id_table == max(id_table))[1])
        
        genesymbol
        
    })
    
    unlist(results)
}