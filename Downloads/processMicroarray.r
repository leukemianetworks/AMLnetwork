###########
##Example##
###########
# pd <- 'PATH_TO_SOURCE_CODE'
# > setwd(pd) # set directory to the folder containing this code
# > source('processMicroarray.r')

# (Starting from raw microarray data with CEL format) 
# wd <- 'PATH_TO_RAW_CEL_DATA'
# > setwd(wd) # Set directory to the folder containing the raw CEL file
# > expr <- readrawdata('GSE30375_RAW.tar','hgu133a.db')
# 						File Name			Library Name



# Load platform database.
# If platform database is not previously downloaded, it will automatically try to download it

#####################
##Loading libraries##
#####################
source("http://bioconductor.org/biocLite.R")
libs <- c('R.utils','tcltk','hash','annotate','simpleaffy','Biobase','GEOquery','affyPLM','limma')
for (lib in libs) {
    # Check if library is installed
    if (eval(parse(text=paste('require(',lib,')',sep='')))) {
        # If yes, load the library and proceed to next
        eval(parse(text=paste('library(',lib,')',sep='')))
        next
    } else {
        # If no, try to install library from CRAN
        eval(parse(text=paste('install.packages("',lib,'")',sep='')))
    }
    # After trying to get library from CRAN, check if the library is installed
    if (eval(parse(text=paste('require(',lib,')',sep='')))) {
        # If yes, load the lirbrary and proceed to next
        eval(parse(text=paste('library(',lib,')',sep='')))
        next
    } else {
        # If no, try to install library from bioconductor
        eval(parse(text=paste('biocLite("',lib,'")',sep='')))
        eval(parse(text=paste('library(',lib,')',sep='')))
    }
    # After tying to get library from both CRAN and biocondutor, if the library is still missing, exit with error message
    if (eval(parse(text=paste('!require(',lib,')',sep='')))) {
        stop('Library could not be found in CRAN nor bioconductor. Please check the supported R version of the library.')
    }
}
message('All required libraries loaded')

#########################
##Function: readrawdata##
#########################
readrawdata <- function(fname, database, files = 'all') {
    # Read raw data (usually in CEL format) and return the normalized data.
    # Example:
    # > pd <- 'PATH_TO_SOURCE_CODE' # pd is the directory where you store this source code
    # > source('processMicroarray.r')
    # > wd <- 'PATH_TO_RAW_DATA' # wd is the directory where you put your raw data
    # > setwd('wd') 
    # > expr <- readrawdata('GSE30375_RAW.tar','hgu133a.db') 
	#						Filename			library name
	# Currently, the name of raw data is default to be GSEXXXXX_RAW.tar, which is also the default raw data name from GEO database.
    
    # Load platform database.
    # If platform database is not previously downloaded, it will automatically try to download it

    options(warn=-1)
    Sys.sleep(0.1)
    message(paste('Loading platform database: ',database,sep=''))
    if (database %in% installed.packages()) {
        eval(parse(text=paste('library(',database,')',sep='')))
    } else {
        message('Unable to find platform database. Trying to download.')
        eval(parse(text=paste('biocLite("',database,'")',sep='')))
        if (database %in% installed.packages()) {
            eval(parse(text=paste('library(',database,')',sep='')))
        } else {
            stop('Unable to find/download platform database. Please check platform name.')
        }
    }
    options(warn=0)
    Sys.sleep(0.1)
    message('Platform database loaded')

    # Load raw microarray data
    unlink('data', recursive=TRUE)
    Sys.sleep(0.1)
    message('Loading raw microarray data')
    gse <- strsplit(fname,'_')[[1]][1]

    # Raw data should be untar to data folder. 
    untar(fname,exdir='data/')
    # Each CEL data should then be compressed in gz format.
    cels <- list.files('data/', pattern = 'gz')
    sapply(paste('data', cels, sep='/'), gunzip)
    file.remove(paste('data', cels, sep='/'))
    # Read data
    cels <- paste('data', list.files('data/', pattern='.CEL'), sep='/')
    if (class(files) == 'numeric' && length(files) == 2) {
        cels <- cels[files[1]:files[2]]
    }
    celData <- ReadAffy(filenames=cels)
    Sys.sleep(0.1)
    message('Raw microarray data loaded')

    # Normalization & Write to file
    # threestep
    # http://www.bioconductor.org/packages/release/bioc/vignettes/affyPLM/inst/doc/ThreeStep.pdf
    Sys.sleep(0.1)
    message('Normalizing data')
    eset <- threestep(celData)
    Sys.sleep(0.1)
    message('Data normalized')
    # Write normalized expression in tab-delimited txt file and compressed using gzip
    mat <- exprs(eset)
    condition <- colnames(mat)
    genes <- rownames(mat)
    mat <- rbind(condition,mat)
    rownames(mat)[1] <- gse
    colnames(mat) <- NULL
    gzf <- gzfile(paste(gse,'_NormalizedFromRaw.txt.gz',sep=''),'w')
    write.table(mat, gzf, sep="\t", row.name=TRUE, col.names=FALSE, quote=FALSE)
    close(gzf)

    # Clean up & Return
    unlink('data', recursive=TRUE)
    expr <- list()
    expr$gse <- gse
    expr$database <- database
    expr$data <- exprs(eset)

	
	annotatedata(expr)
   # return(expr)
}


##########################
##Function: annotatedata##
##########################
annotatedata <- function(expr) {
    # Annotate data and mapping from probes to gene with median value
    # Example:
    # > annotatedata(expr) # expr is the output from readrawdata

    # Load platform database.
    # If platform database is not previously downloaded, it will automatically try to download it
    database <- expr$database
    options(warn=-1)
    Sys.sleep(0.1)
    message(paste('Loading platform database: ',database,sep=''))
    if (database %in% installed.packages()) {
        eval(parse(text=paste('library(',database,')',sep='')))
    } else {
        message('Unable to find platform database. Trying to download.')
        eval(parse(text=paste('biocLite("',database,'")',sep='')))
        if (database %in% installed.packages()) {
            eval(parse(text=paste('library(',database,')',sep='')))
        } else {
            stop('Unable to find/download platform database. Please check platform name.')
        }
    }
    options(warn=0)
    Sys.sleep(0.1)
    message('Platform databse loaded')
	message('annotating data')

    # Map gene symbols
    gene.symbols <- getSYMBOL(rownames(expr$data), expr$database)
    results <- cbind(gene.symbols, expr$data)
    # Remove NA symbols
    results <- results[!isNA(gene.symbols),]
    gene.symbols <- results[,'gene.symbols']
    # Probe to Gene Mapping with Median
    genes <- unique(gene.symbols)
    results.gene <- c()
	#WINDOWS ONLY    pb <- tkProgressBar(title = "progress bar", min = 0, max = length(genes))
	pb <- txtProgressBar(title = "progress bar", min = 0, max = length(genes))
	
    c <- 1
    Sys.sleep(0.1)
    message('Progress Bar')
    for (g in 1:length(genes)) {
        gene <- genes[g]
        result <- results[gene.symbols==gene,]
        if (is.null(nrow(result))) {
            results.gene <- rbind(results.gene,result)
        } else {
            tmp <- c(gene)
            for (i in 2:dimension(result)[2]) {
                tmp <- c(tmp,median(as.numeric(result[,i])))
            }
            results.gene <- rbind(results.gene,tmp)
        }
        Sys.sleep(0.1)
        #WINDOWS ONLY  setTkProgressBar(pb, c, label=paste( round(c/length(genes)*100, 0),"% done"))
		setTxtProgressBar(pb, g)
        c <- c+1
    }
    close(pb)
    Sys.sleep(0.1)
    message('Data annotated')
    message('Writing result')
    #results.gene <- t(results.gene)
    #results.gene <- results.gene[2:dimension(results.gene)[1],]
    colnames(results.gene)[1] <- expr$gse
    rownames(results.gene) <- NULL

    results.gene <- t(results.gene)
    write.table(results.gene, paste('Microarray_',expr$gse,'_Processed.txt',sep=''), sep='\t', row.name=TRUE, col.names=FALSE, quote=FALSE)
    message('Done')
}
