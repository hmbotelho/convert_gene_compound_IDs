# uniprot_to_genename
#
# Converts Uniprot IDs to full text gene names
# 
# Hugo Botelho
# v0.1
# 24 July 2021
#
# Input: 
#	* up: character vector, with Uniprot IDs.
#   * asvector: output as character vector? Othwerwise, the output will be a data.frame.
#   * showProgress: logical, display progress bar in the R console. Disabled when `parallelize = TRUE`
#   * parallelize: logical, use multi-core processing to improve performance.
# Output: character or data.frame, with the input and corresponding gene names
uniprot_to_genename <- function(up, asvector=TRUE, showProgress=TRUE, parallelize=FALSE) {
    
    # Sanity checks
    if(class(up) != "character") stop("Class of 'up' is not 'character'.")
    
    
    if(parallelize == FALSE | parallel::detectCores() < 2){
        # Single-core processing
        
        res <- c()
        
        for(i in seq_along(up)){
            
            if(showProgress){
                txtProgressBar(0, length(up), i, style=3)
            }
            
            tryCatch({
                url        <- paste0("https://www.uniprot.org/uniprot/", up[i])
                htmlsource <- readLines(url, encoding = "UTF-8", warn = FALSE)
                parsed_doc <- XML::htmlParse(htmlsource, encoding = "UTF-8")
                temp       <- XML::xpathSApply(parsed_doc, path = '//*[@id="content-protein"]/h1', XML::xmlValue)
                if(is.null(temp)) temp <- NA
            },
            error = function(e){
                temp <<- NA
            },
            warning = function(w){
                temp <<- NA
            })
            
            res <- c(res, temp)
        }
        
        names(res) <- up
        
    } else {
        # Multi-core processing
        
        # Start cluster
        cluster <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cluster, varlist=c("up"), envir=environment())
        
        # Concatenate objects file
        res <- parLapply(cluster, 1:length(up), function(i){
            
            tryCatch({
                up1        <- up[i]
                url        <- paste0("https://www.uniprot.org/uniprot/", up1)
                htmlsource <- readLines(url, encoding = "UTF-8", warn = FALSE)
                parsed_doc <- XML::htmlParse(htmlsource, encoding = "UTF-8")
                temp       <- XML::xpathSApply(parsed_doc, path = '//*[@id="content-protein"]/h1', XML::xmlValue)
                if(is.null(temp)) temp <- NA
            },
            error = function(e){
                temp <<- NA
            },
            warning = function(w){
                temp <<- NA
            })
            
            c(i, up1, temp)
        })
        
        res <- do.call(rbind, res)
        res <- res[order(as.integer(res[,1])),]
        res <- res[,3]
        names(res) <- up
        
        # Stop cluster
        stopCluster(cluster)
    }
    
    
    # Possible export as data.frame
    if(!asvector){
        res <- data.frame(uniprotswissprot = up,
                          genename         = res,
                          stringsAsFactors = FALSE)
        rownames(res) <- NULL
    }
    
    res
}
