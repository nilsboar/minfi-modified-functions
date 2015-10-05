estimateCellCounts <- function (combLMset, meanPlot = FALSE) {
    # combinedRGset <- combine(rgSet, referenceRGset)
    # newpd <- data.frame(sampleNames = c(sampleNames(rgSet), sampleNames(referenceRGset)),
                        # studyIndex = rep(c("user", "reference"),
                        # times = c(ncol(rgSet), ncol(referenceRGset))),
                        # stringsAsFactors=FALSE)
    # pData(combinedRGset) <- newpd
    # referencePd <- pData(referenceRGset)
    # rm(referenceRGset)
       
    ## Extracts normalized reference data:
    # browser()
    referenceMset <- combLMset[, combLMset$Study == "bloodcell.methy"]
    # pData(referenceMset) <- as(referencePd, "DataFrame") #?  restoring the original pData, already there for us
    lmSet <- combLMset[, combLMset$Study == "dcg.methy"]

    compData <- pickCompProbes(referenceMset)  # fix compData, different object
    coefs <- compData$coefEsts
    
    cat("[estimateCellCounts] Estimating composition.\n")
    data.betas <- estimateBeta(lmSet, returnType = "matrix")
    ref.betas <- estimateBeta(referenceMset, returnType = "matrix")
    counts <- projectCellType(data.betas[rownames(coefs), ], coefs)
    rownames(counts) <- sampleNames(lmSet)  # originally rgSet not lmset

    if (meanPlot) { # original...
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(colMeans(estimateBeta(lmSet)[rownames(coefs), ]), smeans)

        sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
        plot(sampleMeans, pch = 21, bg = sampleColors)
        legend("bottomleft", c("blood", levels(factor(names(smeans)))),
               col = 1:7, pch = 15)
    }
    compTable <- as.matrix(compData$compTable)[,colnames(counts)]
    ######################################
    #  seeking to extract specific CD4 markers
    ######################################
    # look for difference of CD4 column from others:
    
    comp.tbl.ordrd <- compTable[order(compTable[,'CD4T'] - apply(compTable, 1, mean)),]
    max.compTable.cd4 <- comp.tbl.ordrd[470009:470058, ]
    min.compTable.cd4 <- comp.tbl.ordrd[1:50, ]
    # calculate average cell comp across individuals
    avg.cts <- colMeans(counts)  # add to 1.02, fix?
    adj.cts <- sweep(counts, 2, avg.cts)  # subtract the average
    adj.data <- data.betas - compTable %*% t(adj.cts)
    # run through again with adjusted data
    counts.tst <- projectCellType(adj.data[rownames(coefs), ], coefs)
    # now return!
    list(counts = counts, counts.tst = counts.tst, compTable = compData$compTable,  adjData = adj.data, hicd4markers = max.compTable.cd4, locd4markers = min.compTable.cd4)   
}
pickCompProbes <- function (mSet, cellTypes = c("Gran", "CD4T", "CD8T", "Bcell", "Mono", "NK", "Eos", "Neu"), numProbes = 50)  # ? can't get probe list for other three cell types
{
    splitit <- function(x) {
        split(seq(along = x), x)
    }
    browser()
    p <- estimateBeta(mSet, returnType = "matrix")  #  
    pd <- as.data.frame(pData(mSet))
    browser()
    # no use of mSet below here -> we only need p and pd
    if (is.null(cellTypes)) {
        cellTypes <- unique(pd$CellType)
    }
    else {
        if (!all(cellTypes %in% pd$CellType)) 
            stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep, ]
        p <- p[, keep]
    }
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)  # should be ok here  
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[, i]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low", 
        "high", "range")
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0, ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    probeList <- lapply(tstatList, function(x) {
        y <- x[x[, "p.value"] < 1e-08, ]
        yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
        yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
        c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
    })
    trainingProbes <- unlist(probeList)  # trainingProbes looks like empty character vector
    p <- p[trainingProbes, ]  # get error: subscript out of bounds
    pMeans <- rowMeans(p)
    names(pMeans) <- pd$CellType
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), 
        collapse = "+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType - 1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
    out <- list(coefEsts = coefEsts, compTable = compTable, sampleMeans = pMeans)
    return(out)
}

projectCellType  <- function (Y, coefCellType, contrastCellType = NULL, nonnegative = TRUE, 
    lessThanOne = FALSE) 
{
    if (is.null(contrastCellType)) 
        Xmat <- coefCellType
    else Xmat <- coefCellType %*% t(contrastCellType)
    nCol <- dim(Xmat)[2]
    nSubj <- dim(Y)[2]
    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    if (nonnegative) {
        if (lessThanOne) {
            Amat <- cbind(rep(-1, nCol), diag(nCol))
            b0vec <- c(-1, rep(0, nCol))
        }
        else {
            Amat <- diag(nCol)
            b0vec <- rep(0, nCol)
        }
        for (i in 1:nSubj) {
            obs <- which(!is.na(Y[, i]))
            Dmat <- t(Xmat[obs, ]) %*% Xmat[obs, ]
            # browser()
            mixCoef[i, ] <- solve.QP(Dmat, t(Xmat[obs, ]) %*% 
                Y[obs, i], Amat, b0vec)$sol
        }
    }
    else {
        for (i in 1:nSubj) {
            obs <- which(!is.na(Y[, i]))
            Dmat <- t(Xmat[obs, ]) %*% Xmat[obs, ]
            mixCoef[i, ] <- solve(Dmat, t(Xmat[obs, ]) %*% Y[obs, 
                i])
        }
    }
    return(mixCoef)
}

