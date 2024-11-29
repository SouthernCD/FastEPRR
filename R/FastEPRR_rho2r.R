FastEPRR_rho2r <- function(inputFilePath = null, outputFilePath = null, Ne = null) {
    ################# check inputFilePath ###############
    if(is.null(inputFilePath) || !is.character(inputFilePath)) {
        stop("Please input valid 'inputFilePath' (The absolute path of input file)!")
    }
    if(!file.exists(inputFilePath)) {
        stop("Please input valid 'vcfFilePath' (The full path of input file) (ends in [.gz|.vcf])!")
    }
    ################# check inputFilePath ###############

    ################# check outputFilePath ###############
    if(!is.character(outputFilePath)) {
        stop("Please input valid 'srcOutputFilePath' (The absolute path of output file)!")
    }
    if(!file.exists(dirname(outputFilePath))) {
        stop("The 'srcOutputFilePath' is not valid!")
    }
    ################# check outputFilePath ###############

    ################# check outputFilePath ###############
    if(!is.null(Ne) && !is.numeric(Ne)) {
        stop("Please input valid 'Ne' (The effective population size)!")
    }
    if(!is.null(Ne) && Ne <= 0) {
        stop("Please input valid 'Ne' (The effect population size)!")
    }
    ################# check outputFilePath ###############
    
    #data <- read.table(file = inputFilePath, header = TRUE)
    readFile <- file(inputFilePath, "r")
    writeFile <- file(outputFilePath, "w")
    headInfo <- ''
    line <- readLines(readFile, 1)
    n = length(unlist(strsplit(line, split = " ")))
    if(n == 3) {
        headInfo <- c("Start End r(cM/Mb)")
    } else if(n == 5) {
        headInfo <- c("Start End r(cM/Mb) CIL(cM/Mb) CIR(cM/Mb)")
    } else {
        stop("Please input valid file!")
    }
    write(headInfo, file = writeFile)
    line <- readLines(readFile, 1)
    while(length(line) != 0) {
        line = unlist(strsplit(line, split = " "))
        start <- as.numeric(line[1])
        end <- as.numeric(line[2])
        rho <- as.numeric(line[3])
        winLength <- end - start + 1
        r <- rho * 1e8 / (4 * Ne * winLength)
        if(length(line) == 3) {
            write(paste(line[1], line[2], as.character(r)), file = writeFile)
        } else {
            ciL <- as.numeric(line[4])
            ciR <- as.numeric(line[5])
            rL <- ciL * 1e8 / (4 * Ne * winLength)
            rR <- ciR * 1e8 / (4 * Ne * winLength)
            write(paste(line[1], line[2], as.character(r), as.character(rL) ,as.character(rR)), file = writeFile)
            #writeLines(c(start, end ,rho, rL ,rR), con = writeFile)
        }
        line = readLines(readFile, 1)
    }
    close(readFile)
    close(writeFile)
}

