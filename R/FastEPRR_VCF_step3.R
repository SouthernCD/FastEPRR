FastEPRR_VCF_step3 <- function(srcFolderPath = NULL, DXFolderPath = NULL, finalOutputFolderPath = NULL) {

    ################# check srcFilePath ###############
    if(is.null(srcFolderPath) || is.na(srcFolderPath) || !is.character(srcFolderPath)) {
        stop("Please input valid 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
    }
    if(!file.exists(srcFolderPath)) {
        stop("Please input valid 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
    }
    ################# check srcFilePath ###############



    ################# check DXFolderPath ###############
    if(is.null(DXFolderPath) || is.na(DXFolderPath) || !is.character(DXFolderPath)) {
        stop("Please input valid 'DXFolderPath' (The folder path stores files which named by nsam_doulbeton_xton)!")
    }
    if(!file.exists(DXFolderPath)) {
        stop("The 'DXFolderPath' is not valid!")
    }
    ################# check DXFolderPath ###############



    ################# check finalOutputFolderPath ###############
    if(is.null(finalOutputFolderPath) || is.na(finalOutputFolderPath) || !is.character(finalOutputFolderPath)) {
        stop("Please input valid 'finalOutputFolderPath' (The folder path which stores the output files)!")
    }
    if(!file.exists(finalOutputFolderPath)) {
        stop("The 'finalOutputFolderPath' is not exist !")
    }
    ################# check finalOutputFolderPath ###############

    options(scipen = 200)
    allChr <- list.files(path = srcFolderPath, full.names = TRUE)
    if(length(allChr) == 0) {
        stop("Please input valid 'srcFolderPath' (The parent directory of 'srcOutputFilePath' in step1.")
    }

    alldouxData <- list.files(path=DXFolderPath, full.names=TRUE)
    if(length(alldouxData) == 0) {
        stop("Please input valid 'DXFolderPath' (The folder path stores files which named by nsam_doulbeton_xton)!")
    }

    assign("myData", NULL, envir = .GlobalEnv)
    assign("rhoData", NULL, envir = .GlobalEnv)
    assign("hasCI", NULL, envir = .GlobalEnv)

    lapply(allChr, readStep1_step3)

    lapply(alldouxData, readStep2)

    for(i in 1:length(myData)) {
        myData[[i]] <<- cbind(myData[[i]], rhoData[[i]])
    }

    lapply(myData, writeResults, finalOutputFolderPath)

    rm(hasCI, myData, rhoData, envir = .GlobalEnv)
}

readStep1_step3 <- function(fileName) {
    temp <- read.table(file = fileName, header = TRUE, stringsAsFactors = FALSE, numerals = 'no.loss')
    if(nrow(temp) > 0) {
        temp <- temp[, c(1, 10, 11)]
        if(is.null(myData)) {
            myData <<- list(temp)
            names(myData) <<- as.character(temp[1, 1])

            rhoData <<- list(as.data.frame(matrix(nrow = nrow(temp), ncol = 3)))
            names(rhoData) <<- as.character(temp[1, 1])
        } else {
            myData[[length(myData) + 1]] <<- temp
            names(myData)[length(myData)] <<- as.character(temp[1, 1])

            rhoData[[length(rhoData) + 1]] <<- as.data.frame(matrix(nrow = nrow(temp), ncol = 3))
            names(rhoData)[length(rhoData)] <<- as.character(temp[1, 1])
        }
    }
}


readStep2 <- function(fileName) {
    temp <- read.table(file = fileName, header = TRUE, stringsAsFactors =  FALSE)
    # print(fileName)
    if(length(temp[1,]) == 6) {
        hasCI <<- TRUE
    } else {
        hasCI <<- FALSE
    }
    if(nrow(temp) > 0) {
        temp <- as.data.frame(t(temp))
        lapply(temp, allocateDX)
    }

}

allocateDX <- function(dxLine) {
    dxLine <- as.vector(dxLine)
    chrName <- dxLine[1]
    startPos <- as.numeric(dxLine[2]) 
    endPos <- as.numeric(dxLine[3])

    listNum <- grep(paste("^", chrName, "$", sep = ""), names(myData))

    # print(dxLine)
    # print(listNum)
    commons <- which(myData[[listNum]]$startPos == startPos)

    if(length(dxLine) == 4) {
        rhoData[[listNum]][commons, 1] <<- as.numeric(dxLine[4])
    } else {
        rhoData[[listNum]][commons, 1] <<-  as.numeric(dxLine[4])
        rhoData[[listNum]][commons, 2] <<-  as.numeric(dxLine[5])
        rhoData[[listNum]][commons, 3] <<-  as.numeric(dxLine[6])
    }
}

writeResults <- function(singleChr, outputFolderPath) {
    currchr <- singleChr[1,1]
    #hasCI <- length(singleChr[1,]) == 6
    saveFile <- paste(outputFolderPath, "/chr", currchr, sep = "")
    if(file.exists(saveFile)) {
        file.remove(saveFile)
    }
    tSingleChr <- as.data.frame(t(singleChr))
    saveFile <- file(saveFile, "w")
    if(hasCI) {
        writeLines("Start End Rho CIL CIR", con = saveFile)
    } else {
        writeLines("Start End Rho", con = saveFile)
    }
    lapply(tSingleChr, writeSingleChr, saveFile)
    close(saveFile)
}


writeSingleChr <- function(singleItem, fileName) {
    singleItem <- as.vector(singleItem)
    writeLines(paste(as.integer(singleItem[2]), " ", as.integer(singleItem[3]), sep = ""), con = fileName, sep = " ")
    if(is.na(singleItem[4])) {
        stop("Step2 is not accomplished!")
    }
    rho <- as.numeric(singleItem[4])
    rhoL <- as.numeric(singleItem[5])
    rhoR <- as.numeric(singleItem[6])
    if(hasCI) {
        if(rhoL == 0 && rhoR == 0) {
            writeLines(paste(ifelse(rho, rho, 0), " NA", " NA", sep = ""), fileName)
        } else {
            writeLines(paste(ifelse(rho, rho, 0), " ", ifelse(rhoL, rhoL, 0), " ", rhoR, sep = ""), fileName)
        }
    } else {
        writeLines(paste(ifelse(rho, rho, 0)), fileName)
    }
}









