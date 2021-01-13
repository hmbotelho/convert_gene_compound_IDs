# geneid_to_genesymbol
#
# Converts NCBI gene IDs to official gene symbols.
# 
# Hugo Botelho & André Falcão
# v0.2
# 12 January 2021
#
# Input: 
#	* geneids: vector (character or numeric), with NCBI gene IDs
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame
#	* showProgress: logical, print progress to console?
# Output: data frame or character, with the input and corresponding gene symbols
geneid_to_genesymbol <- function(geneids, asvector=TRUE, showProgress = FALSE){
    
    output <- lapply(geneids, function(geneid){
        
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
        
        c(geneid, genesymbol)
        
    })
    
    output <- do.call("rbind", output)
    output <- as.data.frame(output, stringsAsFactors = FALSE)
    
    colnames(output) <- c("entrezgene_id", "hgnc_symbol")
    output <- output[match(geneids, output$entrezgene_id),]
    output$entrezgene_id <- geneids
    class(output$entrezgene_id) <- "integer"
    class(output$hgnc_symbol) <- "character"
    rownames(output) <- 1:nrow(output)
    
    if(asvector){
        output <- output$hgnc_symbol
        names(output) <- geneids
    }
    
    output
}

