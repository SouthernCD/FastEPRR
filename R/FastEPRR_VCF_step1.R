FastEPRR_VCF_step1 <- function(vcfFilePath = NULL, erStart = NULL, erEnd = NULL, winLength = NULL, stepLength = NULL, winDXThreshold = 10, samConsidered = "all", qualThreshold = 20, gapFilePath = NULL, srcOutputFilePath = NULL) {

    #dyn.load("I:/FastEPRR_2.0/R/FastEPRR.dll")

    ################# check vcfFilePath ###############
    if(is.null(vcfFilePath) || !is.character(vcfFilePath)) {
        stop("Please input valid 'vcfFilePath' (The absolute path of input file) (ends in [.gz|.vcf])!")
    }
    if(!file.exists(vcfFilePath) || !grepl(".*[/|\\].*(.gz|.vcf)$", vcfFilePath, ignore.case=TRUE)) {
        stop("Please input valid 'vcfFilePath' (The full path of input file) (ends in [.gz|.vcf])!")
    }
    ################# check vcfFilePath ###############

    ################# check erStart ###############
    options(scipen=200)
    if(!is.null(erStart) && !is.numeric(erStart)) {
        stop("Please input valid 'erStart' (The start position of estimating rho)!")
    }
    if(!is.null(erStart) && erStart <= 0) {
        stop("Please input valid 'erStart' (The start position of estimating rho)!")
    }
    ################# check erStart ###############

    ################# check erEnd ###############
    if(!is.null(erEnd) && !is.numeric(erEnd)) {
        stop("Please input valid 'erEnd' (The start position of estimating rho)!")
    }
    if(!is.null(erEnd) && erEnd <= 0) {
        stop("Please input valid 'erEnd' (The start position of estimating rho)!")
    }
    if(is.numeric(erStart) && is.numeric(erEnd) && (erEnd <= erStart)) {
        stop("The 'erEnd' must be greater than 'erStart'!")
    }
    ################# check erEnd ###############

    ################# check winLength ###############
    if(!is.numeric(winLength)) {
        stop("Please input valid 'winLength' (The length of sliding window)!")
    }
    if(winLength <= 0) {
        stop("Please input valid 'winLength' (The length of sliding window)!")
    }
    ################# check winLength ###############

    ################# check stepLength ###############
    if(!is.null(stepLength)) {
        if(!is.numeric(stepLength)) {
            stop("Please input valid 'stepLength' (The length of sliding step)!")
        }
        if(stepLength <= 0) {
            stop("Please input valid 'stepLength' (The length of sliding step)!")
        }
        if(stepLength > winLength) {
            warning("Because the 'stepLength' is greater than 'winLength', there is no overlapping sliding windows!")
        }
    }
    ################# check stepLength ###############

    ################# check winDXThreshold ###############
    if(!is.numeric(winDXThreshold)) {
        stop("Please input valid 'winDXThreshold' (The minimum number of doubleton and xton of each window)!")
    }
    if(winDXThreshold < 0) {
        stop("Please input valid 'winDXThreshold' (The minimum number of doubleton and xton of each window)!")
    }
    if(winDXThreshold >= 0) {
        assign("winSNPThre", winDXThreshold, envir = .GlobalEnv)
    }
    ################# check winDXThreshold ###############

    ################# check samConsidered ###############
    specifiedIndiNames <- NULL
    specifiedIndiChrs <- NULL
    if(samConsidered == "all") {
        specifiedIndiChrs <- c(1)
    } else if(samConsidered == "all[0:1]") {
        specifiedIndiChrs <- c(2)
    } else if(samConsidered == "all[1:0]") {
        specifiedIndiChrs <- c(3)
    } else {
        allIndis <- FALSE
        if(!is.character(samConsidered)) {
            stop("Please input valid 'samConsidered' (The considered samples)!")
        }
        if(!grepl(";", samConsidered)) {
            stop("Please input valid 'samConsidered' (The considered samples)!")
        }

        samVals <- trim(unlist(strsplit(samConsidered, ";")))

        if(any(is.null(samVals)) || any(is.na(samVals)) || any(nchar(samVals) == 0)) {
            stop("Please input valid 'samConsidered' (The considered samples)!")
        }

        specifiedIndiNames <- vector("character", length = length(samVals))

        if(grepl(":", samConsidered) || grepl("\\[", samConsidered) || grepl("\\]", samConsidered)) {
            commaCount <- length(grep(":", samVals))
            leftBCount <- length(grep("\\[", samVals))
            rightBCount <- length(grep("\\]", samVals))
            if((commaCount != leftBCount) || (leftBCount != rightBCount) || (commaCount != rightBCount)) {
                stop("Please input valid 'samConsidered' (The considered samples)!")
            }

            specifiedIndiChrs <- vector("numeric", length = length(samVals))
            specifiedIndex <- grepl("\\[0:1\\]|\\[1:0\\]", samVals)
            specifiedIndiNames[!specifiedIndex] <- samVals[!specifiedIndex]
            specifiedIndiChrs[!specifiedIndex] <- 1

            posChr <- regexpr("\\[", samVals)
            specifiedIndiNames[specifiedIndex] = substr(samVals[specifiedIndex], 1, posChr[posChr > 0] - 1)
            specifiedIndiChrs[grepl("\\[0:1\\]", samVals)] <- 2
            specifiedIndiChrs[grepl("\\[1:0\\]", samVals)] <- 3

            if(any(specifiedIndiChrs == 0)) {
                stop("Please input valid 'samConsidered' (The considered samples)!")
            }

        } else {
            specifiedIndiNames <- samVals
            specifiedIndiChrs <- rep(1, length(samVals))
        }
    }
    ################# check samConsidered ###############


    ################# check qualThreshold ###############
    if(!is.numeric(qualThreshold)) {
        stop("Please input valid 'qualThreshold' (VCF quality threshold)!")
    }
    ################# check qualThreshold ###############

    ################# check gapFilePath ###############
    if(!is.null(gapFilePath) && (!is.character(gapFilePath))) {
        stop("Please input valid absolute path of 'gapFilePath' (gap file path) or let it be NULL!")
    }
    if(!is.null(gapFilePath) && !file.exists(gapFilePath)) {
        stop("The 'gapFilePath' is not valid!")
    }
    ################# check gapFilePath ###############

    ################# check srcOutputFilePath ###############
    if(!is.character(srcOutputFilePath)) {
        stop("Please input valid 'srcOutputFilePath' (The absolute path of output file)!")
    }
    if(!file.exists(dirname(srcOutputFilePath))) {
        stop("The 'srcOutputFilePath' is not valid!")
    }
    if(file.exists(srcOutputFilePath)) {
        file.remove(srcOutputFilePath)
    }
    ################# check srcOutputFilePath ###############

    srcOutputFilePath <- file(srcOutputFilePath, "w")
    assign("FastEPRRSavePaths", srcOutputFilePath, envir = .GlobalEnv)
    paras <- list(vcfFP = vcfFilePath, erS = erStart, erE = erEnd, winL = winLength, stepL = stepLength, speIndiNames = specifiedIndiNames, speIndiChrs = specifiedIndiChrs, qual = qualThreshold, gapFP = gapFilePath, output = srcOutputFilePath)
    headInfo <- c("chr nsam wdou wxton wH wHetero wavsk2 wavr2 misInfo startPos endPos")
    writeLines(headInfo, con = FastEPRRSavePaths)
    # print(paras)
    vcfEsRho(paras)
    close(srcOutputFilePath)
    rm(FastEPRRSavePaths, winSNPThre, envir = .GlobalEnv)
}

winRhoVCF <- function(currChr, nsam, seqs, leftP, rightP) {
    # print(seqs)

    mfsAndH <- getWinMFSAndH(nsam, t(seqs), nrow(seqs), "VCF")

    if(!is.null(mfsAndH)) {

        fourSSObs <- mfsAndH$fourSSInfo
        mfs <- mfsAndH$mfsInfo
        tempRow <-  NULL

        if(mfsAndH$MDataOrnot) { # have missing data
            tempRow <- c(currChr, nsam, mfs[1], mfs[2], fourSSObs[1], fourSSObs[2], fourSSObs[3], fourSSObs[4], paste(c(mfsAndH$samSize, mfsAndH$MDataIDs), collapse=","), leftP, rightP)
        } else {
            tempRow <- c(currChr, nsam, mfs[1], mfs[2], fourSSObs[1], fourSSObs[2], fourSSObs[3], fourSSObs[4], NA, leftP, rightP)
        }

        writeLines(tempRow, con = FastEPRRSavePaths, sep = " ")
        writeLines("", con = FastEPRRSavePaths)
        rm(mfsAndH)
    }


}
