
vcfEsRho <- function(paras) {
    cat("Parsing file\n\n", paras$vcfFP, "\n\n")
    erS <- paras$erS
    erE <- paras$erE
    if(is.numeric(erS) && is.numeric(erE)) {
        readVCF(paras, erS, erE)
    } else if(is.numeric(erS) && !is.numeric(erE)) {
        readVCF(paras, erS, Inf)
    } else if(!is.numeric(erS) && is.numeric(erE)) {
        readVCF(paras, -Inf, erE)
    } else {
        readVCF(paras, -Inf, Inf)
    }
}

readVCF <- function(paras, rStart, rEnd) {
    gapFile <- NULL
    if(!is.null(paras$gapFP)) {
        tryCatch({
            gapFile <- read.table(paras$gapFP, header = TRUE, colClasses = "numeric", comment.char = "")
            if(ncol(gapFile) > 2) {
                stop("Invalid gap file format!")
            }
        }, error = function(err) {
            stop("'gapFilePath' parameter is not String!")
        }, warning = function(warn) {
            stop("'gapFilePath' parameter is not String!")
        }, finally = {
            
        })
    }
    
    specifiedIndiNames <- paras$speIndiNames
    specifiedIndiChrs <- paras$speIndiChrs
    qualThreshold <- paras$qual
    winLength <- paras$winL
    stepLength <- NULL
    if(is.numeric(paras$stepL)) {
        stepLength <- paras$stepL
    } else {
        stepLength <- winLength
    }
    samConsiderd <- 0
    totalIndiNum <- 0
    usedSamIndex <- NULL
    vcfFile <- file(paras$vcfFP, open = "r")
    CHROMIndex <- 0
    POSIndex <- 0
    REFIndex <- 0
    ALTIndex <- 0
    QUALIndex <- 0
    FORMATIndex <- 0
    leftP <- rStart
    rightP <- leftP + winLength - 1
    dataCopy <- NULL
    posCopy <- NULL
    checkPhase <- TRUE
    while(length(line <- readLines(vcfFile, n = 3000, warn = FALSE)) > 0) {
        #print(line)
        if(grepl("/", line[1])) {
            stop("This VCF file has unphased genotype data!\n")
        }
        headInfo <- unlist(strsplit(line[grepl("#CHR", line)], " |\t"))
        if(length(headInfo) != 0) {
            CHROMIndex <- grep("#CHROM", headInfo)
            POSIndex <- grep("POS", headInfo)
            REFIndex <- grep("REF", headInfo)
            ALTIndex <- grep("ALT", headInfo)
            QUALIndex <- grep("QUAL", headInfo)
            FORMATIndex <- grep("FORMAT", headInfo)
            # print(CHROMIndex)
            # print(POSIndex)
            # print(REFIndex)
            # print(ALTIndex)
            # print(QUALIndex)
            # print(FORMATIndex)
            if((CHROMIndex != 1) ||POSIndex != 2 || REFIndex != 4 || ALTIndex != 5 || QUALIndex != 6|| FORMATIndex != 9) {
                stop("This VCF file has a wrong format, please check!\n")
            }
            totalIndiNum <- length(headInfo) - FORMATIndex
            
            if(is.null(specifiedIndiNames)) {
                ## Index of used samples
                usedSamIndex <- 10:(totalIndiNum * 2 + FORMATIndex)
                if(specifiedIndiChrs == 1) {
                    
                } else if(specifiedIndiChrs == 2) {
                    usedSamIndex <- usedSamIndex[usedSamIndex %% 2 == 1]
                } else {
                    usedSamIndex <- usedSamIndex[usedSamIndex %% 2 != 1]
                }
            } else {
                usedSamIndex <- vector("numeric", length = 0)
                for(i in 1:length(specifiedIndiNames)) {
                    s <- grep(specifiedIndiNames[i], headInfo)
                    if(length(s) == 0) {
                        stop("There are some individuals not in VCF file !\n")
                    } else if(length(s) > 1) {
                        stop("There are same individuals in VCF file !\n")
                    }
                    if(specifiedIndiChrs[i] == 1) {
                        usedSamIndex <- c(usedSamIndex, (s - 10) * 2 + 10, (s - 10) * 2 + 11)
                    } else if(specifiedIndiChrs[i] == 2) {
                        usedSamIndex <- c(usedSamIndex, (s - 10) * 2 + 11)
                    } else {
                        usedSamIndex <- c(usedSamIndex, (s - 10) * 2 + 10)
                    }
                }
            }
            samConsiderd <- length(usedSamIndex) 
            if(samConsiderd < 6) {
                stop("The sample size required must be greater than or equal to six!")
            }
            if(samConsiderd <= floor(log(samConsiderd))) {
                stop(cat("The minimum value of winLength for sample size ", samConsiderd, " is ", (floor(log(samConsiderd)) + 1) ,"!", sep=""))
            }
            headInfo <- NULL
        }
        if(checkPhase) {
            checkPhase <- FALSE
            if(grepl("/", line[1])) {
                stop("This VCF file has unphased genotype data !\n")
            }
        }
        matrixData <- matrix(unlist(strsplit(line[!grepl("#", line)], " |\\||\t")), ncol = totalIndiNum * 2 + FORMATIndex, byrow = TRUE)
        ## Get rid of positions < rStart
        startPosIndex <- which(as.numeric(matrixData[, POSIndex]) > rStart)
        if(length(startPosIndex) == 0) {
            next()
        }
        ## Get rid of positions > rEnd
        endPosIndex <- which(as.numeric(matrixData[, POSIndex]) < rEnd)
        if(length(endPosIndex) == 0) {
            break()
        }
        ## Get rid of indels and low quality data
        readIndex <- intersect(snpIndex(matrixData[, REFIndex], matrixData[, ALTIndex]), qualIndex(matrixData[, QUALIndex], qualThreshold))
        rangeIndex <- intersect(startPosIndex, endPosIndex)
        readIndex <- intersect(readIndex, rangeIndex)
        ## Postions of used data
        usedPos <- as.numeric(matrixData[readIndex, POSIndex])
        ## Chr of used data
        usedChr <- matrixData[1, CHROMIndex][1]
        
        ## Filter data
        usedData <- matrixData[readIndex, usedSamIndex, drop = FALSE]
        if(!is.null(posCopy)) {
            usedPos <- c(posCopy, usedPos)
            posCopy <- NULL
        }
        if(!is.null(dataCopy)) {
            usedData <- rbind(dataCopy, usedData)
            dataCopy <- NULL
        }
        rm(matrixData)
        ## Check data
        # lapply(usedData, checkAlleles)
        
        if(leftP == -Inf) {
            leftP <- usedPos[1]
            rightP <- leftP + winLength - 1
        }
        # print("#######")
        # print(usedPos)
        # print("#######")
        if(length(usedPos) == 0) {
            next()
        }
        while(rightP < usedPos[length(usedPos)]) {
            if(!insecGap(leftP, rightP, gapFile)) {
                calIndex <- intersect(which(usedPos >= leftP), which(usedPos <= rightP))
                if(length(calIndex) != 0) {
                    if(length(calIndex) == 1) {
                        winRhoVCF(usedChr, samConsiderd, t(usedData[calIndex,]), leftP, rightP)
                        # print(calIndex)
                        # print(leftP)
                        # print(rightP)
                    } else {
                        winRhoVCF(usedChr, samConsiderd, usedData[calIndex,], leftP, rightP)
                        # print(calIndex)
                        # print(leftP)
                        # print(rightP)
                    }
                }
            }
            leftP <- leftP + stepLength
            rightP <- rightP + stepLength
            if(rightP > rEnd) {
                rightP <- rEnd
                break()
            }
        }
        calIndex <- which(usedPos >= leftP)
        if(length(calIndex) == 0) {
            posCopy <- NULL
            dataCopy <- NULL
        } else {
            posCopy <- usedPos[calIndex]
            dataCopy <- usedData[calIndex,]
        }
    }
    close(vcfFile)
    if(!is.null(dataCopy)) {
        if(length(dataCopy) == samConsiderd) {
            winRhoVCF(usedChr, samConsiderd, t(dataCopy), leftP, usedPos[length(usedPos)])
            # print(leftP)
            # print(rightP)
        } else {
            winRhoVCF(usedChr, samConsiderd, dataCopy, leftP, usedPos[length(usedPos)])
            # print(leftP)
            # print(rightP)
        }
    }
}

snpIndex <- function(refAllele, altAllele) {
    ref <- which(nchar(refAllele) > 1)
    altPoint <- which(grepl("\\.", altAllele))
    altLTO <- which(nchar(altAllele) > 1)
    altSNP <- vector("numeric", length = 0)
    alt <- NULL
    if(length(altLTO) != 0) {
        for(i in 1:length(altLTO)) {
            if(grepl(",", altAllele[altLTO[i]]) && all(grepRaw(",", altAllele[altLTO[i]], all = TRUE) %% 2 == 0)) {
                altSNP <- c(altSNP, i)
            }
        }
        altLTO <- setdiff(altLTO, altSNP)
        alt <- union(altPoint, altLTO)
    } else {
        alt <- altPoint
    }
    indel <- union(ref, alt)
    snp <- setdiff(1:length(refAllele), indel)
    return(snp)
}

qualIndex <- function(quality, threshold) {
    pointIndex <- which(quality == ".")
    numberIndex <- which(as.numeric(quality[quality != "."]) >= threshold)
    return(union(pointIndex, numberIndex))
}

insecGap <- function(leftP, rightP, gapData) {
    condition1 <- gapData[,1] < rightP
    condition2 <- leftP <= gapData[,2]
    if(length(intersect(which(condition1 == TRUE), which((condition2 == TRUE))))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

checkAlleles <- function(alleles) {
    if(alleles == "0" || alleles == "1" || alleles == "2" || alleles == "\\.") {
        
    } else {
        stop("There are some characters falling outside the range of '0', '1', '2' and '.'!\n")
    }
}















