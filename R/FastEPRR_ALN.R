FastEPRR_ALN <- function(alnFilePath = NULL, fileFormat = NULL, erStart = NULL, erEnd = NULL, winLength = NULL, stepLength = NULL, getCI = TRUE, replicateNum = 100, demoParameter = NULL, trainingSet1 = NULL, trainingSet2 = NULL, outputFilePath = NULL) {

    #dyn.load("I:/FastEPRR_2.0/R/FastEPRR.dll")

    ################# check vcfFilePath ###############
    if(is.null(alnFilePath) || !is.character(alnFilePath)) {
        stop("Please input valid 'alnFilePath' (The absolute path of input file) (ends in [.fas|.aln|.phy])!")
    }
    if(!file.exists(alnFilePath) || !grepl(".*[/|\\].*(.fas|.aln|.phy)$", alnFilePath, ignore.case=TRUE)) {
        stop("Please input valid 'alnFilePath' (The full path of input file) (ends in [.fas|.aln|.phy])!")
    }
    ################# check vcfFilePath ###############


    ################# check fileFormat ###############
    if(is.null(fileFormat) || !is.numeric(fileFormat)) {
        stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")
    }
    if(fileFormat == 1) {
        if(!grepl(".*[/|\\].*(.fas)$", alnFilePath)) {
            stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")
        } else {
            fileFormat <- "fasta"
        }
    } else if(fileFormat == 2) {
        if(!grepl(".*[/|\\].*(.aln)$", alnFilePath)) {
            stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")
        } else {
            fileFormat <- "clustal"
        }
    } else if(fileFormat == 3) {
        if(!grepl(".*[/|\\].*(.phy)$", alnFilePath)) {
            stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")
        } else{
            fileFormat <- "phylip"
        }
    } else {
        stop("Please input int value 'format' (format number) of alignment file  (1:fasta; 2:clustal; 3:phylip)!")
    }
    ################# check fileFormat ###############

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
    if(!is.null(winLength) && !is.numeric(winLength)) {
        stop("Please input valid 'winLength' (The length of sliding window)!")
    }
    if(is.numeric(winLength) && winLength <= 0) {
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

    ################# check replicateNum ###############
    if(!is.numeric(replicateNum)) {
        stop("Please input valid 'replicateNum'!")
    }

    assign("FastEPRRrepN", replicateNum, envir = .GlobalEnv)
    ################# check replicateNum ###############

    ################# check getCI ###############
    if(!is.logical(getCI)) {
        stop("Please input valid 'getCI' (whether calculate CI for each window)!")
    }
    assign("needCI", getCI, envir = .GlobalEnv)
    ################# check getCI ###############

    ################# outputFilePath demoParameter ###############
    if(is.null(outputFilePath) || !is.character(outputFilePath)) {
        stop("Please input valid character value 'outputFilePath' (The absolute path of output file)!")
    }
    if(!file.exists(dirname(outputFilePath))) {
        stop("The 'outputFilePath' is not valid!")
    }
    if(file.exists(outputFilePath)) {
        file.remove(outputFilePath)
    }
    ################# outputFilePath demoParameter ###############

    ################# check trainingSet1 ###############
    if(!is.null(trainingSet1)) {
        tryCatch({
            if(!is.character(trainingSet1)) {
                stop("Please input valid 'training set' (the range rho of training set)!")
            }
            num <- unlist(strsplit(trainingSet1, ";"))
            num <- as.numeric(num)
            assign("rhoTrainingSet1", num, envir = .GlobalEnv)
        }, error = function(err) {
            stop("Please input valid 'training set' (the range rho of training set)!")
        }, warning = function(warn) {
            stop("Please input valid 'training set' (the range rho of training set)!")
        }, finally = {

        })
    } else {
        assign("rhoTrainingSet1", c(0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 70.0, 110.0, 170.0), envir = .GlobalEnv)
    }
    ################# check trainingSet1 ###############

    ################# check trainingSet2 ###############
    if(!is.null(trainingSet2)) {
        tryCatch({
            if(!is.character(trainingSet2)) {
                stop("Please input valid 'training set' (the range rho of training set)!")
            }
            num <- unlist(strsplit(trainingSet2, ";"))
            num <- as.numeric(num)
            assign("rhoTrainingSet2", num, envir = .GlobalEnv)
        }, error = function(err) {
            stop("Please input valid 'training set' (the range rho of training set)!")
        }, warning = function(warn) {
            stop("Please input valid 'training set' (the range rho of training set)!")
        }, finally = {

        })

        if(max(rhoTrainingSet1) >= min(rhoTrainingSet2)) {
            stop("Rhos in training set 2 must be larger than rhos in training set 1!")
        }
    } else {
        assign("rhoTrainingSet2", c(180.0, 190.0, 200.0, 220.0, 250.0, 300.0, 350.0), envir = .GlobalEnv)
    }
    ################# check trainingSet2 ###############

    ################# check demoParameter ###############
    if(!is.null(demoParameter) && (!is.character(demoParameter))) {
        stop("Please input valid character value 'demoParameter' (The population demographic parameters)!")
    }
    if(!is.null(demoParameter)) {
        demoInfo <- unlist(strsplit(demoParameter, " "))
        if((length(which(demoInfo == "-f")) > 0) || (length(which(demoInfo == "-seeds")) > 0) || (length(which(demoInfo == "-t")) > 0) ||
           (length(which(demoInfo == "-s")) > 0) || (length(which(demoInfo == "-T")) > 0) || (length(which(demoInfo == "-L")) > 0) ||
           (length(which(demoInfo == "-p")) > 0) || (length(which(demoInfo == "-r")) > 0) || (length(which(demoInfo == "-c")) > 0)) {
            stop("'-f -seeds -t -s -T -L -p -r -c' parameters are not allowed!")
        }
        if((length(which(demoInfo == "f")) > 0) || (length(which(demoInfo == "seeds")) > 0) || (length(which(demoInfo == "t")) > 0) ||
           (length(which(demoInfo == "s")) > 0) || (length(which(demoInfo == "T")) > 0) || (length(which(demoInfo == "L")) > 0) ||
           (length(which(demoInfo == "p")) > 0) || (length(which(demoInfo == "r")) > 0) || (length(which(demoInfo == "c")) > 0)) {
            stop("'-f -seeds -t -s -T -L -p -r -c' parameters are not allowed!")
        }
        assign("FastEPRRdemoParas", demoParameter, envir = .GlobalEnv)
    } else {
        assign("FastEPRRdemoParas", NULL, envir = .GlobalEnv)
    }
    ################# check demoParameter ###############


    installPkes("seqinr")

    paras <- list(fp = alnFilePath, fm = fileFormat, es = erStart, ee = erEnd, wl = winLength, sl = stepLength)

    assign("winSNPThre", NULL, envir = .GlobalEnv)

    assign("FastEPRRDataInfos", data.frame(nsam = integer(), wdou = integer(), wxton = integer(), wH = integer(), wHetero = numeric(), wavsk2 = numeric(), wavr2 = numeric(),
        misInfo = character(), startPos = integer(), endPos = integer(), rho = numeric(), rhoL = numeric(), rhoR = numeric(), flag = integer(), stringsAsFactors = FALSE), envir = .GlobalEnv)

    alignEsRho(paras)
    startBoosting(outputFilePath)
    rm(FastEPRRDataInfos, FastEPRRdemoParas, FastEPRRrepN, needCI, rhoTrainingSet1, rhoTrainingSet2, winSNPThre, envir = .GlobalEnv)
}

startBoosting <- function(outputFilePath) {
    installPkes("mboost")

    outputFilePath <- file(outputFilePath, "w")
    countWin <- 1
    cat("Total windows: ", nrow(FastEPRRDataInfos), "\n")

    if(needCI) {
        writeLines("Start End Rho CIL CIR", con = outputFilePath)
    } else {
        writeLines("Start End Rho", con = outputFilePath)
    }

    while(countWin <= nrow(FastEPRRDataInfos)) {

        cat("window size:", countWin, "\n")
        currWin <- FastEPRRDataInfos[countWin, ]

        writeLines(paste(currWin$startPos, " ", currWin$endPos, sep = ""), sep = " ", con = outputFilePath)


        currDou <- currWin$wdou
        currXto <- currWin$wxton
        currH <- currWin$wH

        douchkMis <- TRUE

        if(is.na(currWin$misInfo)) {
            currNsam <- currWin$nsam
            missDaIDs <- integer(0)
        } else { # have missing data
            missAndsam <- as.integer(unlist(strsplit(currWin$misInfo, ",")))
            currNsam <- missAndsam[1]
            missDaIDs <- missAndsam[-1]
            if(length(missDaIDs) == 0) {
                douchkMis <- FALSE
            }
        }

        if(is.na(currWin$flag)) {
            if((length(missDaIDs)==0) & douchkMis) {
                sameConfigALN(currNsam, currDou, currXto)
                currWin <- FastEPRRDataInfos[countWin, ]
                if(needCI) {
                    if(currWin$rhoR == 0.00 || currWin$rhoR == "0.00") {
                        writeLines(paste(currWin$rho, "  NA", " NA", sep=""), con = outputFilePath)
                    } else {
                        writeLines(paste(currWin$rho, " ", currWin$rhoL, " ", currWin$rhoR, sep = ""), con = outputFilePath)
                        # cat("Rho:", rhoOutputFor(as.numeric(currWin$rho)), " CIL:", rhoOutputFor(as.numeric(currWin$rhoL)), " CIR:", rhoOutputFor(as.numeric(currWin$rhoR)), "\n", sep="", file=outputFilePath, append=TRUE)
                    }
                } else {
                    writeLines(paste(currWin$rho), con = outputFilePath)
                }
            } else { # have missing data
                mfsAndSS <- list(samSize = currNsam, mfsInfo = c(currDou, currXto), fourSSInfo = c(currWin$wH, currWin$wHetero, currWin$wavsk2, currWin$wavr2), MDataIDs = missDaIDs)
                erhos <- missEstiRho(mfsAndSS)
                if(needCI) {
                    if(erhos[2] == 0 & erhos[2] == "0.00") {
                        writeLines(paste(erhos[3], " NA", " NA", sep=""), con = outputFilePath)
                    } else {
                        writeLines(paste(erhos[3], " ", erhos[1], " ", erhos[2], sep = ""), con = outputFilePath)
                        # cat("Rho:", rhoOutputFor(as.numeric(currWin$rho)), " CIL:", rhoOutputFor(as.numeric(currWin$rhoL)), " CIR:", rhoOutputFor(as.numeric(currWin$rhoR)), "\n", sep="", file=outputFilePath, append=TRUE)
                    }
                } else {
                    writeLines(paste(erhos[3]), con = outputFilePath)
                }
            }
        } else {
            if(needCI) {
                if(currWin$rhoR == 0.00 & currWin$rhoR == "0.00") {
                    writeLines(paste(currWin$rho, " NA", " NA", sep=""), con = outputFilePath)
                } else {
                    writeLines(paste(currWin$rho, " ", currWin$rhoL, " ", currWin$rhoR, sep = ""), con = outputFilePath)
                    # cat("Rho:", rhoOutputFor(as.numeric(currWin$rho)), " CIL:", rhoOutputFor(as.numeric(currWin$rhoL)), " CIR:", rhoOutputFor(as.numeric(currWin$rhoR)), "\n", sep="", file=outputFilePath, append=TRUE)
                }
            } else {
                writeLines(paste(currWin$rho), con = outputFilePath)
            }

        }
        countWin <- countWin + 1
    }
    close(outputFilePath)
}

sameConfigALN <- function(nsam, doubleton, xton) {

    rhoIndex <- 11
    rhoLIndex <- 12
    rhoRIndex <- 13
    flagIndex <- 14

    missDaIDs <- integer(0)
    wZero <- getZeroLowHigh(nsam, doubleton, xton, missDaIDs)
    zeroLow <- wZero[1]
    zeroHigh <- wZero[2]

    bmlAndH <- firstModel(nsam, doubleton, xton, missDaIDs)

    thrH <- bmlAndH$thrH
    firBml <- bmlAndH$firMod
    secBml <- NULL

    commons <- which((FastEPRRDataInfos$wdou == doubleton) & (FastEPRRDataInfos$wxton == xton) & is.na(FastEPRRDataInfos$misInfo))
    for(i in commons) {
        currWin <- FastEPRRDataInfos[i, ]
        currH <- currWin$wH
        if(currH <= zeroLow) {
            FastEPRRDataInfos[i, rhoIndex] <<- 0.00
            FastEPRRDataInfos[i, rhoLIndex] <<- 0.00
            FastEPRRDataInfos[i, rhoRIndex] <<- 0.00
            FastEPRRDataInfos[i, flagIndex] <<- 1
        } else if(currH <= thrH) {
            mfsAndSS <- list(samSize = currWin$nsam, mfsInfo = c(doubleton, xton), fourSSInfo = c(currWin$wH, currWin$wHetero, currWin$wavsk2, currWin$wavr2), MDataIDs = integer(0))
            returnRho <- alphaCorrect(mfsAndSS, firBml)

            if(needCI) {
                CIul <- getCI(mfsAndSS, firBml, returnRho, rhoTrainingSet1[length(rhoTrainingSet1)])
                upCI <- max(CIul)
                if((currH >= zeroLow) && (currH <= zeroHigh)) {
                    downCI <- 0.0
                } else {
                    downCI <- min(CIul)
                }
                if(downCI < 0.0) {
                    downCI <- 0
                }
                FastEPRRDataInfos[i, rhoIndex] <<- returnRho
                FastEPRRDataInfos[i, rhoLIndex] <<- downCI
                FastEPRRDataInfos[i, rhoRIndex] <<- upCI
                FastEPRRDataInfos[i, flagIndex] <<- 1
            } else {
                FastEPRRDataInfos[i, rhoIndex] <<- returnRho
                FastEPRRDataInfos[i, rhoLIndex] <<- 0.00
                FastEPRRDataInfos[i, rhoRIndex] <<- 0.00
                FastEPRRDataInfos[i, flagIndex] <<- 1
            }

            rm(mfsAndSS)
        } else {
            mfsAndSS <- list(samSize = currWin$nsam, mfsInfo = c(doubleton, xton), fourSSInfo = c(currWin$wH, currWin$wHetero, currWin$wavsk2, currWin$wavr2), MDataIDs = integer(0))

            if(is.null(secBml)) {
                secBml <- secModel(nsam, doubleton, xton, missDaIDs, bmlAndH$firTD)
            }

            returnRho <- alphaCorrect(mfsAndSS, secBml)

            if(needCI) {
                CIul <- getCI(mfsAndSS, secBml, returnRho, rhoTrainingSet2[length(rhoTrainingSet2)])
                upCI <- max(CIul)
                if((currH >= zeroLow) && (currH <= zeroHigh)) {
                    downCI <- 0.0
                } else {
                    downCI <- min(CIul)
                }
                if(downCI < 0.0) {
                    downCI <- 0
                }
                FastEPRRDataInfos[i, rhoIndex] <<- returnRho
                FastEPRRDataInfos[i, rhoLIndex] <<- downCI
                FastEPRRDataInfos[i, rhoRIndex] <<- upCI
                FastEPRRDataInfos[i, flagIndex] <<- 1
            } else {
                FastEPRRDataInfos[i, rhoIndex] <<- returnRho
                FastEPRRDataInfos[i, rhoLIndex] <<- 0.00
                FastEPRRDataInfos[i, rhoRIndex] <<- 0.00
                FastEPRRDataInfos[i, flagIndex] <<- 1
            }
            rm(mfsAndSS)
        }
    }
}



alnWinRho <- function(nsam, seqs, leftP, rightP) {

    mfsAndH <- getWinMFSAndH(nsam, seqs, ncol(seqs), "ALN")

    if(!is.null(mfsAndH)) {
        fourSSObs <- mfsAndH$fourSSInfo
        mfs <- mfsAndH$mfsInfo
        tempRow <- c(nsam, mfs[1], mfs[2], fourSSObs[1], fourSSObs[2], fourSSObs[3], fourSSObs[4])
        if(mfsAndH$MDataOrnot) { # have missing data
            tempRow <- c(tempRow, paste(c(mfsAndH$samSize, mfsAndH$MDataIDs), collapse=","), leftP, rightP, NA, NA, NA, NA)
        } else {
            tempRow <- c(tempRow, NA, leftP, rightP, NA, NA, NA, NA)
        }
        dInfosRow <- nrow(FastEPRRDataInfos) + 1
        FastEPRRDataInfos[dInfosRow, ] <<- tempRow
    }

}
