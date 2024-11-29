FastEPRR_VAR <- function(varFilePath = NULL, winLength = NULL, stepLength = NULL, gapFilePath = NULL, varOutputFilePath = NULL) {

    ################# check varFilePath ###############
    if(is.null(varFilePath) || is.na(varFilePath) || !is.character(varFilePath)) {
        stop("Please input 'varFilePath' (The absolute path of input file)!")
    }
    if(!file.exists(varFilePath)) {
        stop("Please input 'varFilePath' (The absolute path of input file)!")
    }
    ################# check varFilePath ###############


    ################# check winLength ###############
    if(!is.numeric(winLength)) {
        stop("Please input valid 'winLength' (The length of sliding window)!")
    }
    if(winLength <= 0) {
        stop("Please input valid 'winLength' (The length of sliding window)!")
    }
    ################# check winLength ###############


    ################# check stepLength ###############
    if(!is.numeric(stepLength)) {
        stop("Please input valid 'stepLength' (The length of sliding step)!")
    }
    if(winLength / stepLength != 2) {
        stop("'winLength' is not the double of 'stepLength'!")
    }
    ################# check stepLength ###############


    ################# check gapFilePath ###############
    if(!is.null(gapFilePath) && (!is.character(gapFilePath))) {
        stop("Please input valid absolute path of 'gapFilePath' (gap file path) or let it be NULL!")
    }
    if(!is.null(gapFilePath) && !file.exists(gapFilePath)) {
        stop("The 'gapFilePath' is not valid!")
    }
    ################# check gapFilePath ###############


    ################# check varOutputFilePath ###############
    if(is.null(varOutputFilePath) || is.na(varOutputFilePath) || !is.character(varOutputFilePath)) {
        stop("Please input valid character value 'varOutputFilePath' (The absolute path of output file)!")
    }
    if(!file.exists(dirname(varOutputFilePath))) {
        stop("The 'varOutputFilePath' is not valid!")
    }
    if(file.exists(varOutputFilePath)) {
        file.remove(varOutputFilePath)
    }
    assign("saveFile", file(varOutputFilePath, "w"), envir = .GlobalEnv)
    ################# check varOutputFilePath ###############

    writeLines("Start End Rho", con = saveFile)
    options(scipen=200)

    assign("posAndRhoList", list(), envir = .GlobalEnv)

    readFileVar(varFilePath, winLength, stepLength, gapFilePath)

    lapply(posAndRhoList, obtainRho, stepLength)

    close(saveFile)
    rm(forwPos, posAndRhoList, saveFile, segCC, envir = .GlobalEnv)
}

readFileVar <- function(filePath, winlength, stepLength, gapFilePath) {

    gapFile <- NULL
    if(!is.null(gapFilePath)) {
        tryCatch({
            gapFile <- read.table(gapFilePath, header = TRUE, colClasses = "numeric", comment.char = "")
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

    assign("forwPos", NULL, envir = .GlobalEnv)
    proWinLen <- NULL

    inFile <- file(filePath, open="r")
    inBreak <- FALSE
    assign("segCC", 1, envir = .GlobalEnv)

    startIndex <- 1
    endIndex <- 2
    rhoIndex <- 3

    lines <- readLines(inFile, n = -1, warn = FALSE)
    numOfline <- length(unlist(strsplit(lines[1], " ")))
    lines <- lines[-1]
    matrixData <- matrix(unlist(strsplit(lines, " ")), ncol = numOfline, byrow = TRUE)
    if(is.null(proWinLen)) {
        proWinLen <- as.numeric(matrixData[1, endIndex]) - as.numeric(matrixData[1, startIndex]) + 1
        if(proWinLen != winlength) {
            stop("Your input 'winLength' is not equal to the window length in the input file !")
        }
    }
    matrixData <- as.data.frame(t(matrixData[,c(startIndex, endIndex, rhoIndex)]), stringsAsFactors = FALSE)
    lapply(matrixData, readInLine, stepLength, gapFile)

    close(inFile)
    rm(list = ls())
}

readInLine <- function(aline, stepLength, gapData) {
    startIndex <- 1
    endIndex <- 2
    rhoIndex <- 3

    startPos <- as.numeric(aline[startIndex])
    endPos <- as.numeric(aline[endIndex])
    rho <- as.numeric(aline[rhoIndex])
    inBreak <- FALSE

    if(is.null(forwPos)) {
        posAndRhoList[[segCC]] <<- list()
        posAndRhoList[[segCC]][[1]] <<- list(lp = startPos, rp = endPos, r = rho)
    } else {
        if(forwPos + stepLength == startPos) {
            tlen <- length(posAndRhoList[[segCC]]) + 1
            posAndRhoList[[segCC]][[tlen]] <<- list(lp = startPos, rp = endPos, r = rho)
        } else if(forwPos + stepLength > startPos) {
            stop("The positions are not correct, please check it!")
        } else {
            tempPos <- forwPos + stepLength
            while(tempPos < startPos) {
                if(insecGap(tempPos, tempPos + 2*stepLength - 1 ,gapData)) {
                    segCC <<- segCC + 1
                    inBreak <- TRUE
                    break
                } else {
                    tlen <- length(posAndRhoList[[segCC]]) + 1
                    posAndRhoList[[segCC]][[tlen]] <<- list(lp = tempPos, rp = tempPos + 2*stepLength - 1, r = NA)
                }
                tempPos <- tempPos + stepLength
            }

            if(inBreak) {
                posAndRhoList[[segCC]] <<- list()
            }

            tlen <- length(posAndRhoList[[segCC]]) + 1
            posAndRhoList[[segCC]][[tlen]] <<- list(lp = startPos, rp = endPos, r = rho)
        }
    }
    forwPos <<- startPos
}

obtainRho <- function(rhoList, stepLength) {
    rhoMatrix <- matrix(unlist(rhoList), ncol = 3, byrow = TRUE)
    naIndex <- which(is.na(rhoMatrix[,3]))

    if(length(naIndex) != 0) {
        #print(length(naIndex))
        #print(naIndex)
        #print(rhoMatrix)
        if(is.na(rhoMatrix[length(rhoMatrix)])) {
            rNa <- is.na(rhoMatrix[,3])
            count <- length(rNa)
            while(count > 0) {
                if(!rNa[count]) {
                    rhoMatrix <- rhoMatrix[1:count,, drop = FALSE]
                    naIndex <- which(is.na(rhoMatrix[,3]))
                    break
                }
                count <- count - 1
            }
        }
        #print(naIndex)
        if(length(naIndex) != 0) {

            startRho <- as.numeric(rhoMatrix[naIndex[1] - 1, 3])
            endRho <- 0
            startP <- naIndex[1]
            if(length(naIndex) != 1) {

                for(i in 1:(length(naIndex) - 1)) {
                    if(naIndex[i] + 1 != naIndex[i + 1]) {
                        endRho <- as.numeric(rhoMatrix[naIndex[i] + 1, 3])
                        rhoMatrix[startP:naIndex[i], 3] <- (startRho + endRho) / 2
                        startRho <- as.numeric(rhoMatrix[naIndex[i + 1] - 1, 3])
                        startP <- naIndex[i + 1]
                    }
            }
            }
            endRho <- as.numeric(rhoMatrix[naIndex[length(naIndex)] + 1, 3])
            rhoMatrix[startP:naIndex[length(naIndex)], 3] <- (startRho + endRho) / 2

        }
    }

    len <- nrow(rhoMatrix)
    comIndex <- 1
    x1 <- 0
    x2 <- 0
    x3 <- 0
    x4 <- 0

    c2 <- 0
    c3 <- 0
    c4 <- 0

    L <- 0
    U <- 0
    A <- 0
    B <- 0
    C <- 0
    rhosList <- list()
    verySmall <- 0.00001
    f3 <- function(x, a, b, c, d) {
        a*x^3 + b*x^2 + c*x +d
    }
    while((comIndex+2) <= len) {
        # print(rhoMatrix[comIndex,])
        rho1 <- as.numeric(rhoMatrix[comIndex, 3])
        rho2 <- as.numeric(rhoMatrix[comIndex + 1, 3])
        rho3 <- as.numeric(rhoMatrix[comIndex + 2, 3])

        if(rho2 >= (rho1 + rho3)) {
            x1 <- 0
            x4 <- 0
            tmp <- (rho2 - rho1 - rho3)/3
            x2 <- tmp + rho1
            x3 <- tmp + rho3
        } else {
            c2 <- rho1
            c3 <- rho2 - rho1
            c4 <- rho3 - rho2 + rho1
            L <- max(0, -c3)
            U <- min(c2, c4)
            A <- c3 - c2 - c4
            B <- c2*c4 - c3*c4 - c3*c2
            C <- c2*c3*c4
            rootSet <- c(L, U)
            tryCatch({
                realRoot <- uniroot(f3, lower = L + verySmall, upper = U - verySmall, a = 4, b = 3*A, c = 2*B, d = C, tol = 0.00001)$root
                rootSet <- c(rootSet, realRoot)
            }, error = function(err) {

            }, warning = function(warn) {

            }, finally = {

            })
            yValue <- NULL
            for(i in rootSet) {
                yValue <- c(yValue, fx(i, c2, c3, c4))
            }
            x1 <- rootSet[yValue == max(yValue)][1]
            x2 <- c2 - x1
            x3 <- c3 + x1
            x4 <- c4 - x1
        }
        if(x1 >= 0 && x2 >= 0 && x3 >= 0 && x4 >= 0) {
            x1lp <- as.numeric(rhoMatrix[comIndex, 1])
            rhosList <- storeRho(rhosList, x1lp, x1lp + stepLength - 1, x1)

            x2lp <- as.numeric(rhoMatrix[comIndex + 1, 1])
            rhosList <- storeRho(rhosList, x2lp, x2lp + stepLength - 1, x2)

            x3lp <- as.numeric(rhoMatrix[comIndex + 2, 1])
            rhosList <- storeRho(rhosList, x3lp, x3lp + stepLength - 1, x3)

            x4lp <- x3lp + stepLength
            rhosList <- storeRho(rhosList, x4lp, x4lp + stepLength - 1, x4)
        }

        comIndex <- comIndex + 1
    }
    lapply(rhosList, saveFineRho)

}

saveFineRho <- function(singleList) {
    writeLines(paste(singleList$lp, " ", singleList$rp, " ", mean(singleList$r), sep = ""), con = saveFile)
}

storeRho <- function(rhosList, leftP, rightP, rho) {
    len <- length(rhosList)
    if(len == 0) {
        rhosList[[1]] <- list(lp = leftP, rp = rightP, r = rho)
    } else {
        findIndex <- -1
        for(i in 1:len) {
            if(rhosList[[i]]$lp == leftP && rhosList[[i]]$rp == rightP) {
                findIndex <- i
                break
            }
        }

        if(findIndex != -1) { # has find
            tmpRhos <- c(rhosList[[findIndex]]$r, rho)
            rhosList[[findIndex]] <- list(lp = leftP, rp = rightP, r = tmpRhos)
        } else { # not find
            rhosList[[len + 1]] <- list(lp = leftP, rp = rightP, r = rho)
        }
    }

    return(rhosList)
}

fx <- function(x, c2, c3, c4) {
    return(x*(c2-x)*(c3+x)*(c4-x))
}









