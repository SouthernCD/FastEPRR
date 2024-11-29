FastEPRR_VCF_step2 <- function(srcFolderPath = NULL, jobNumber = 1, currJob = 1, demoParameter = NULL, replicateNum = 100, getCI = TRUE, trainingSet1 = NULL, trainingSet2 = NULL, DXOutputFolderPath = NULL) {
    #dyn.load("I:/FastEPRR_2.0/test/FastEPRR.dll")

    ################# check srcFilePath ###############
    if(is.null(srcFolderPath) || is.na(srcFolderPath) || !is.character(srcFolderPath)) {
        stop("Please input valid 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1)!")
    }
    if(!file.exists(srcFolderPath)) {
        stop("Please input valid 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1)!")
    }
    ################# check srcFilePath ###############

    ################# check jobNumber  ###############
    tryCatch({
        if(is.null(jobNumber ) || is.na(jobNumber ) || !is.numeric(jobNumber )) {
            stop("Please input valid 'jobNumber ' (The number of jobs)!")
        }
    }, error=function(err) {
        stop("Please input valid 'jobNumber ' (The number of jobs)!")
    }, warning=function(warn) {
        stop("Please input valid 'jobNumber ' (The number of jobs)!")
    },finally={

    })
    ################# check jobNumber  ###############

    ################# check currJob ###############
    tryCatch({
        if(is.null(currJob) || is.na(currJob) || !is.numeric(currJob)) {
            stop("Please input valid 'currJob' (The current job number)!")
        }
    }, error=function(err) {
        stop("Please input valid 'currJob' (The current job number)!")
    }, warning=function(warn) {
        stop("Please input valid 'currJob' (The current job number)!")
    },finally={

    })
    if(currJob > jobNumber ) {
        stop("currJob must be less than or equal to jobNumber !")
    }
    ################# check currJob ###############

    ################# check demoParameter ###############
    if(!is.null(demoParameter) && (!is.character(demoParameter))) {
        stop("Please input valid 'demoParameter' (The population demographic parameters)!")
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

    ################# check DXOutputFolderPath ###############
    if(is.null(DXOutputFolderPath) || is.na(DXOutputFolderPath) || !is.character(DXOutputFolderPath)) {
        stop("Please input valid 'DXOutputFolderPath' (The folder path stores files which named by DX_doubleton_xton)!")
    }
    if(!file.exists(DXOutputFolderPath)) {
        stop("The 'DXOutputFolderPath' is not valid!")
    }
    assign("outputPath", DXOutputFolderPath, envir = .GlobalEnv)
    ################# check DXOutputFolderPath ###############

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
    options(scipen = 200)

    allChrs <- list.files(path=srcFolderPath, full.names=TRUE)
    if(length(allChrs) == 0) {
        stop("Please input valid character value 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
    }

    assign("allSrcData", NULL, envir = .GlobalEnv)

    lapply(allChrs, readStep1)

    wuniMFS <- allSrcData[, c(2,3,4)]
    colnames(wuniMFS) <- c("nsam", "douton", "xton")
    uniqueMFS <- unique(wuniMFS)
    uniqRow <- nrow(uniqueMFS)

    rm(wuniMFS)

    if(uniqRow < jobNumber) {
        stop(paste("The maximum value of jobNumber  is ", uniqRow, sep=""))
    }

    cat("Total job: ", jobNumber , "\n", "Current job: ", currJob, "\n", sep="")

    eachFileDXs <- floor(uniqRow/jobNumber)

    currDXs <- NULL

    leftWin <- (currJob - 1) * eachFileDXs + 1
    rightWin <- NULL
    if(currJob == jobNumber) {
        rightWin <- uniqRow
    } else {
        rightWin <- currJob * eachFileDXs
    }

    currDXs <- as.data.frame(t(uniqueMFS[leftWin:rightWin,]))

    installPkes("mboost")

    totalDXsIn <- ncol(currDXs)

    assign("saveFile", NULL, envir = .GlobalEnv)
    assign("currentDX", 1, envir = .GlobalEnv)


    lapply(currDXs, dealDX, totalDXsIn)

    allSrcData <<- NULL
    rm(allSrcData, currentDX, FastEPRRdemoParas, FastEPRRrepN, needCI, outputPath, rhoTrainingSet1, rhoTrainingSet2, saveFile, envir = .GlobalEnv)
}

dealDX <- function(currDX, totalDX) {
    cat("Total DXs in this job: ", totalDX, ", current DX: ", currentDX, "\n", sep = "")
    currentDX <<- currentDX + 1

    currNsam <- currDX[1]
    currDou <- currDX[2]
    currXto <- currDX[3]

    saveFile <- paste(outputPath, "/DX", "_", currDou, "_", currXto, sep = "")

    if(file.exists(saveFile)) {
        file.remove(saveFile)
    }

    headInfo <- NULL
    if(needCI) {
        headInfo <- c("chrName startPos endPos rho rhoL rhoR")
    } else {
        headInfo <- c("chrName startPos endPos rho")
    }

    currFile <- file(saveFile, "w")
    writeLines(headInfo, con = currFile)
    cat("Curr dton: ", currDou, " xton: ", currXto, "\n", sep="")

    sameConfig(currNsam, currDou, currXto, currFile)

    close(currFile)
}


readStep1 <- function(fileName) {
    temp <- NULL
    tryCatch({
        temp <- read.table(file = fileName, header = TRUE, stringsAsFactors = FALSE)
    }, error=function(err) {
        stop(paste("File", fileName, "is not valid!"))
    }, warning=function(warn) {
        stop(paste("File", fileName, "is not valid!"))
    },finally={

    })
    allSrcData <<- rbind(allSrcData, temp)
}

sameConfig <- function(nsam, doubleton, xton, resultFile) {

    # Index in results of step1 (input of step2)
    chrIndex <- 1
    nsamIndex <- 2
    wdouIndex <- 3
    wxtonIndex <- 4
    wHIndex <- 5
    wHeteroIndex <- 6
    wavskIndex <- 7
    wavr2Index <- 8
    misInfoIndex <- 9
    startPosIndex <- 10
    endPosIndex <- 11

    missDaIDs <- integer(0)
    assign("secBml", NULL, envir = .GlobalEnv)

    commonWins <- allSrcData[(allSrcData[,wdouIndex] == doubleton) & (allSrcData[,wxtonIndex] == xton) & is.na(allSrcData[,misInfoIndex]),]
    commonWins <- as.data.frame(t(commonWins), stringsAsFactor=FALSE) # for using lapply
    if(length(commonWins) != 0) {
        wZero <- getZeroLowHigh(nsam, doubleton, xton, missDaIDs)
        bmlAndH <- firstModel(nsam, doubleton, xton, missDaIDs)

        lapply(commonWins, estimateRho, wZero, bmlAndH, resultFile)
    }
    # for missing data
    commonWins <- allSrcData[(allSrcData[,wdouIndex] == doubleton) & (allSrcData[,wxtonIndex] == xton) & !is.na(allSrcData[,misInfoIndex]),]
    commonWins <- as.data.frame(t(commonWins), stringsAsFactor=FALSE)
    if(length(commonWins) != 0) {
        lapply(commonWins, missingRho, resultFile)
    }

    rm(secBml, envir = .GlobalEnv)
}


estimateRho <- function(aWin, zeroRange, firM, resultFile) {

    # Index in results of step1 (input of step2)
    chrName <- aWin[1]
    nsamIndex <- 1
    wdouIndex <- 2
    wxtonIndex <- 3
    wHIndex <- 4
    wHeteroIndex <- 5
    wavskIndex <- 6
    wavr2Index <- 7
    misInfoIndex <- 8
    startPosIndex <- 9
    endPosIndex <- 10

    zeroLow <- zeroRange[1]
    zeroHigh <- zeroRange[2]

    thrH <- firM$thrH
    firBml <- firM$firMod

    aWin <- aWin[seq(2, length(aWin))]
    if(!is.numeric(chrName)) {
        suppressWarnings(aWin <- as.numeric(levels(aWin))[aWin])
    }

    currH <- aWin[wHIndex]
    if(currH <= zeroLow) {
        if(needCI) {
            writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], "0.00", "0.00", "0.00", sep = " "), con = resultFile)
            # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], "0.00", "0.00", "0.00", "\n", file = resultFile, append = TRUE)
        } else {
            writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], "0.00", sep = " "), con = resultFile)
            # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], "0.00", "\n", file = resultFile, append = TRUE)
        }
    } else if(currH <= thrH) {
        mfsAndSS <- list(samSize = aWin[nsamIndex], mfsInfo = c(aWin[wdouIndex], aWin[wxtonIndex]), fourSSInfo = c(aWin[wHIndex], aWin[wHeteroIndex], aWin[wavskIndex], aWin[wavr2Index]), MDataIDs = integer(0))

        returnRho <- alphaCorrect(mfsAndSS, firBml)

        if(needCI) {
            CIul <- getCI(mfsAndSS, firBml, returnRho, rhoTrainingSet1[length(rhoTrainingSet1)])
            upCI <- max(CIul)
            downCI <- 0.0
            if((currH >= zeroLow) && (currH <= zeroHigh)) {
                downCI <- 0.0
            } else {
                downCI <- min(CIul)
            }
            if(downCI < 0.0) {
                downCI <- 0
            }

            writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], returnRho, downCI, upCI, sep = " "), con = resultFile)
            # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], rhoOutputFor(returnRho), rhoOutputFor(downCI), rhoOutputFor(upCI), "\n", file = resultFile, append = TRUE)

        } else {
            writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], returnRho, sep = " "), con = resultFile)
            # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], rhoOutputFor(returnRho), "\n", file = resultFile, append = TRUE)
        }

        rm(mfsAndSS)

    } else {
        mfsAndSS <- list(samSize = aWin[nsamIndex], mfsInfo=c(aWin[wdouIndex], aWin[wxtonIndex]), fourSSInfo=c(aWin[wHIndex], aWin[wHeteroIndex], aWin[wavskIndex], aWin[wavr2Index]), MDataIDs=integer(0))

        if(is.null(secBml)) {
            secBml <<- secModel(aWin[nsamIndex], aWin[wdouIndex], aWin[wxtonIndex], integer(0), firM$firTD)
        }

        returnRho <- alphaCorrect(mfsAndSS, secBml)

        if(needCI) {
            CIul <- getCI(mfsAndSS, secBml, returnRho, rhoTrainingSet2[length(rhoTrainingSet2)])
            upCI <- max(CIul)
            downCI <- 0.0
            if((currH >= zeroLow) && (currH <= zeroHigh)) {
                downCI <- 0.0
            } else {
                downCI <- min(CIul)
            }
            if(downCI < 0.0) {
                downCI <- 0
            }

            writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], returnRho, downCI, upCI, sep = " "), con = resultFile)
            # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], rhoOutputFor(returnRho), rhoOutputFor(downCI), rhoOutputFor(upCI), "\n", file = resultFile, append = TRUE)

        } else {
            writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], returnRho, sep = " "), con = resultFile)
            # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], rhoOutputFor(returnRho), "\n", file = resultFile, append = TRUE)
        }

        rm(mfsAndSS)
    }
}

missingRho <- function(aWin, resultFile) {
    # Index in results of step1 (input of step2)
    chrName <- aWin[1]
    nsamIndex <- 1
    wdouIndex <- 2
    wxtonIndex <- 3
    wHIndex <- 4
    wHeteroIndex <- 5
    wavskIndex <- 6
    wavr2Index <- 7
    misInfoIndex <- 8
    startPosIndex <- 9
    endPosIndex <- 10

    aWin <- aWin[seq(2, length(aWin))]

    aWin <- as.character(aWin)
    missAndsam <- as.numeric(unlist(strsplit(aWin[misInfoIndex], ",")))
    if(!is.numeric(chrName)) {
        suppressWarnings(aWin <- as.numeric(aWin))
    }

    currNsam <- missAndsam[1]
    missDaIDs <- missAndsam[-1]
    mfsAndSS <- list(samSize = aWin[nsamIndex], mfsInfo = c(aWin[wdouIndex], aWin[wxtonIndex]), fourSSInfo = c(aWin[wHIndex], aWin[wHeteroIndex], aWin[wavskIndex], aWin[wavr2Index]), MDataIDs = missDaIDs)

    erhos <- missEstiRho(mfsAndSS)

    if(needCI) {
        writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], erhos[1], erhos[2], erhos[3], sep = " "), con = resultFile)
        # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], rhoOutputFor(erhos[1]), rhoOutputFor(erhos[2]), rhoOutputFor(erhos[3]), "\n", file = resultFile, append = TRUE)
    } else {
        writeLines(paste(chrName, aWin[startPosIndex], aWin[endPosIndex], erhos[1], sep = " "), con = resultFile)
        # cat(aWin[chrIndex], aWin[startPosIndex], aWin[endPosIndex], rhoOutputFor(erhos[1]), "\n", file = resultFile, append = TRUE)

    }

}



















