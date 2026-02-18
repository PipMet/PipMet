#' Normalize data
#'
#' This function normalize data.
#' @export
#' @param myDir Path to the directory of work.
#' @param pslist List of spectra.
#' @param metadata A matrix or data.frame with metadata information about samples. Include, at least 'sample' and 'file' columns with name of sample and its path, respectively. More information can be added in new columns, such as 'group', 'class', 'biorep' and 'tecrep'.
#' @param anIC A 'xsAnnotate' CAMERA object with grouped spectra.
#' @param group Character. Name of metadata column to group the samples for normalization. Default to NULL.
#' @param pre_anno A table with annotations for spectra. or path to .csv file.
#' @param removeCompounds Logical. If TRUE, the user may choose from a pop-up identified compounds to remove from the quantification table. If NULL, the user will be asked. Default to 'NULL'.
#' @param mergeCompounds Logical. If TRUE, the user may choose from a pop-up what compounds identified should be representated as one with intensities summed. Exemple: derivatizations derivatives. If NULL, the user will be asked. Default to 'NULL'.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param derivatization Character. Kind of derivatization the samples were prepared with. Supported are 'Trimethylsilyl' and 'None'. If NULL, the user will be asked. Default to 'NULL'.
#' @return A matrix of all spectra, with their annotation (if available), most intense peak m/z and its intensities in every sample.
#' @importFrom grDevices dev.off pdf png tiff boxplot.stats
#' @importFrom graphics boxplot grid legend par text
#' @importFrom svDialogs dlg_message dlg_list
#' @importFrom methods as new
#' @importFrom stats cor sd t.test var
#' @importFrom utils choose.dir menu read.csv select.list write.csv write.table
#' @importFrom tools file_path_as_absolute
#' @import NormalyzerDE
#' @import ggplot2
#' @examples
#' \donttest{
#' \dontrun{
#' pre_anno <- read.csv(system.file('extdata', 'pre_anno.csv', package = 'PipMet'))
#' load(system.file('extdata', 'pslist.RData', package = 'PipMet'))
#' load(system.file('extdata', 'anIC.RData', package = 'PipMet'))
#' load(system.file('extdata', 'metadata.RData', package = 'PipMet'))
#' normalized <- normalize_data(
#'   anIC,
#'   pslist,
#'   metadata,
#'   myDir = '~/',
#'   pre_anno = system.file('extdata', 'pre_anno.csv', package = 'PipMet'),
#'   derivatization = 'Trimethylsilyl',
#'   mergeCompounds = FALSE,
#'   removeCompounds = FALSE,
#'   group = 'group'
#' )
#' }
#' }
normalize_data <- function(anIC, pslist, metadata, myDir, pre_anno, pic_extension = c(".tiff",
    ".png"), derivatization = NULL, mergeCompounds = NULL, removeCompounds = NULL,
    group = NULL) {

    if (!is.null(pre_anno)) {
        if (is_path(pre_anno)) {
            r <- read.csv(pre_anno, sep = ",", na.string = c("NA", ""))
            if (!"annotation" %in% colnames(r)) {
                r <- read.csv(pre_anno, sep = ";", na.string = c("NA", ""), dec = ",")
            }
        } else {
            r <- pre_anno
        }
    } else {
        r <- read.csv("pre_anno.csv", sep = ",", na.string = c("NA", ""))
        if (!"annotation" %in% colnames(r)) {
            r <- read.csv("pre_anno.csv", sep = ";", na.string = c("NA", ""),
                dec = ",")
        }
    }
    r[is.na(r)] <- ""
    pre_anno <- r

    # check if representative ions are ok
    okay <- 2
    while (okay == 2) {
        setwd(myDir)
        quant <- matrix(nrow = nrow(pre_anno), ncol = 5 + (2 * nrow(metadata)))
        colnames(quant) <- c("id", "Fragment Ion (m/z Quant)", "Compound Name",
            "Chemical Formula", "Metabolic Class", metadata$sample, metadata$sample)
        quant[, "id"] <- as.numeric(pre_anno[, "id"])
        quant[, "Compound Name"] <- pre_anno[, "annotation"]

        # set to folder
        if (!dir.exists("Statistics") == TRUE) {
            dir.create("Statistics")
        }
        setwd("Statistics")

        # ask about derivatizations
        if (is.null(derivatization)) {
            derivatization <- dlg_list(c("Trimethylsilyl", "None"), multiple = FALSE,
                title = "Sort of derivatization: ")$res
        }


        for (i in seq_len(nrow(quant))) {
            temp <- anIC@pspectra[[as.integer(quant[i, 1])]]
            if (derivatization == "Trimethylsilyl" && 73 %in% round(pslist[[i]]@spectrum[,
                1])) {
                # if derivatization='Trimethylsilyl', ion m/z 73 is present
                x <- sort(pslist[[i]]@spectrum[-which(round(pslist[[i]]@spectrum[,
                  1]) == 73), 2], decreasing = TRUE)  # remove m/z 73 from possibilities of representative ion
            } else {
                x <- sort(pslist[[i]]@spectrum[, 2], decreasing = TRUE)
            }
            y <- rbind(anIC@groupInfo[temp, ])
            z <- rbind(y[which(y == x[1], arr.ind = TRUE)[, 1], ], y[which(y ==
                x[2], arr.ind = TRUE)[, 1], ])
            quant[i, 2] <- paste0(z[1, 1], ", ", z[2, 1])
            quant[i, 6:(5 + nrow(metadata))] <- z[1, which(colnames(z) == "X1"):which(colnames(z) ==
                paste0("X", nrow(metadata)))]
            quant[i, (6 + nrow(metadata)):(5 + 2 * nrow(metadata))] <- z[2, which(colnames(z) ==
                "X1"):which(colnames(z) == paste0("X", nrow(metadata)))]
        }

        # merge any identified compound (ex: double identified)
        if (is.null(mergeCompounds)) {
            a <- dlg_message("Merge compounds?", type = "yesno")$res
        } else {
            if (mergeCompounds == TRUE) {
                a <- "yes"
            }
            if (mergeCompounds == FALSE) {
                a <- "no"
            }
        }
        if (a == "yes") {
            okay <- "yes"
            while (okay == "yes") {
                z <- dlg_list(sort(quant[which(!is.na(quant[, "Compound Name"]) &
                  !quant[, "Compound Name"] == ""), "Compound Name"]), multiple = TRUE,
                  title = "Merge")$res
                if (length(z) >= 2) {
                  # the user must choose at least two compounds to merge
                  y <- vector()
                  for (i in seq_along(unique(z))) {
                    y <- rbind(y, quant[which(quant[, "Compound Name"] == unique(z)[i]),
                      ])
                  }
                  x <- rbind(y[which(y[, "id"] == sort(as.integer(y[, "id"]))[1]),
                    ])
                  for (i in 6:ncol(quant)) {
                    x[1, i] <- sum(unlist(lapply(y[, i], as.numeric)))
                  }
                  quant[which(quant[, "id"] == sort(as.integer(y[, "id"]))[1]),
                    ] <- x[1, ]
                  for (i in 2:length(z)) {
                    # remove the row of the others identified compounds
                    # after summing to the first
                    quant <- quant[-which(quant[, "id"] == sort(as.integer(y[,
                      "id"]))[i]), ]
                  }
                  okay <- dlg_message("Repeat?", type = "yesno")$res
                } else {
                  okay <- "no"
                }
            }
        }

        # remove any compound identified?
        if (is.null(removeCompounds)) {
            a <- dlg_message("Remove compounds?", type = "yesno")$res
        } else {
            if (removeCompounds == TRUE) {
                a <- "yes"
            }
            if (removeCompounds == FALSE) {
                a <- "no"
            }
        }
        if (a == "yes") {
            okay <- "yes"
            while (okay == "yes") {
                z <- dlg_list(sort(quant[which(!is.na(quant[, "Compound Name"]) &
                  !quant[, "Compound Name"] == ""), "Compound Name"]), multiple = TRUE,
                  title = "Remove")$res
                quant <- quant[-which(quant[, "Compound Name"] == unique(z)[i]),
                  ]
                okay <- dlg_message("Repeat?", type = "yesno")$res
            }
        }

        # write non normalized data into .csv file
        write.csv(quant, "NotNormalized_quantification.csv", row.names = FALSE,
            na = "")

        if (nrow(quant) > 1 & nrow(metadata) > 1) {
            # get the quantification data
            mat <- quant[, 6:ncol(quant)]
            mat <- apply(mat, c(1, 2), FUN = as.numeric)

            # Normalyzer: evaluation and picking of normalization method
            # for the spectra intensities (only most intense peak)
            x <- as.data.frame(mat)
            x[x == 0] <- NA
            x[x > 0 & x < 1] <- 0
            write.table(x[, seq_len(nrow(metadata))], "data.tsv", sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
            write.table(metadata, "design.tsv", sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE)

            # set up for NormalyzerDE
            designFp <- file_path_as_absolute("design.tsv")
            dataFp <- file_path_as_absolute("data.tsv")

            if (is.null(group)) {
                # ask condition to compare from the metadata table
                group <- dlg_list(colnames(metadata), multiple = FALSE, title = "Condition to group from:")$res
            }
            if (nrow(mat) > 50) {
                normalyzer(jobName = "Normalyzer_results", designPath = designFp,
                  dataPath = dataFp, outputDir = myDir, sampleColName = "sample",
                  groupColName = group, requireReplicates = FALSE)
            } else {
                try(normalyzer(jobName = "Normalyzer_results", designPath = designFp,
                  dataPath = dataFp, outputDir = myDir, sampleColName = "sample",
                  groupColName = group, requireReplicates = FALSE, sampleAbundThres = nrow(mat)))
                try(dev.off(), silent = TRUE)  # if the features are less than a threshold, the normalyzer function will crash in the middle of report generating
            }

            # Pick method and apply
            fill <- list.files(paste0(myDir, "/Normalyzer_results"), full.names = TRUE,
                pattern = "-normalized.txt", recursive = TRUE)
            norms <- sub("-normalized.txt", replacement = "", fixed = TRUE, x = basename(fill))

            # ask method for normalization.
            bestNormMat <- menu(norms, graphics = TRUE, title = "Choose the best normalization")

            mat[mat == 0] <- NA
            mat[mat > 0 & mat < 1] <- 0
            if (norms[bestNormMat] == "CycLoess") {
                mat <- performCyclicLoessNormalization(mat)
            }
            if (norms[bestNormMat] == "GI") {
                mat <- globalIntensityNormalization(mat)
            }
            if (norms[bestNormMat] == "log2") {
                mat <- log2(mat)
            }
            if (norms[bestNormMat] == "mean") {
                mat <- meanNormalization(mat)
            }
            if (norms[bestNormMat] == "median") {
                mat <- medianNormalization(mat)
            }
            if (norms[bestNormMat] == "Quantile") {
                mat <- performQuantileNormalization(mat)
            }
            if (norms[bestNormMat] == "RLR") {
                mat <- performGlobalRLRNormalization(mat)
            }
            if (norms[bestNormMat] == "VSN") {
                mat <- performVSNNormalization(mat)
            }
            n <- as.data.frame(cbind(quant[, seq_len(5)], mat))

            # normalize by another criteria (ex: mass, number of cells)
            g <- dlg_message("Normalize by mass/numer of cells/etc?", "yesno")$res
            if (g == "yes") {
                write.csv(n, "Normalized_quantification.csv", row.names = FALSE,
                  na = "")
                dlg_message("Normalize the \"Normalized_quantification.csv\" file as you wish and press \"ok\" to continue.",
                  "ok")$res
                n <- read.csv("Normalized_quantification.csv", na = "", check.names = FALSE)
            }


            # Confirm if the chosen ions are indeed from the same spectrum
            # if everything is ok, use only first ion for the rest of
            # statistics
            y <- data.frame()
            for (ii in 6:(5 + nrow(metadata))) {
                for (i in seq_len(nrow(n))) {
                  y[i, (ii - 5)] <- as.double(n[i, ii])/as.double(n[i, (ii +
                    nrow(metadata))])
                }
            }
            y[, nrow(metadata) + 1] <- apply(y[, seq_len(nrow(metadata))], MARGIN = 1,
                FUN = var, na.rm = TRUE)
            y[, nrow(metadata) + 2] <- apply(y[, seq_len(nrow(metadata))], MARGIN = 1,
                FUN = sd, na.rm = TRUE)
            rownames(y) <- n$id
            colnames(y) <- c(metadata$sample, "Variance", "sd")
            a <- ggplot(y, aes(x = y[, nrow(metadata) + 1])) + geom_histogram() +
                ggtitle("Variance Distribution") + xlab("Variance") + ylab("Frequency") +
                theme(rect = element_rect(fill = "transparent"))  # histograma da variância. O ideal é que tenha pouca variância alta
            b <- ggplot(y, aes(x = y[, nrow(metadata) + 2])) + geom_histogram() +
                ggtitle("Standard Deviation") + xlab("Standard Deviation") +
                ylab("Frequency") + theme(rect = element_rect(fill = "transparent"))  # histograma do desvio padrão.

            # plot histogram of variance and standard deviation and
            # boxplots png
            if (".png" %in% pic_extension) {
                png("variance_dp.png", units = "cm", width = 16, height = 16,
                  res = 900, bg = "NA")
                gridExtra::grid.arrange(a, b, ncol = 1)
                dev.off()
                png("boxplot.png", units = "cm", width = 16, height = 16, res = 900,
                  bg = "NA")
                par(mfrow = c(2, 1))
                boxplot(y$Variance, main = "Variance", horizontal = TRUE)
                boxplot(y$sd, main = "Standard Deviation", horizontal = TRUE)
                dev.off()
            }

            # tiff
            if (".tiff" %in% pic_extension) {
                tiff("variance_dp.tiff", units = "cm", width = 16, height = 16,
                  res = 900, bg = "NA")
                gridExtra::grid.arrange(a, b, ncol = 1)
                dev.off()
                tiff("boxplot.tiff", units = "cm", width = 16, height = 16, res = 900,
                  bg = "NA")
                par(mfrow = c(2, 1))
                boxplot(y$Variance, main = "Variance", horizontal = TRUE)
                boxplot(y$sd, main = "Standard Deviation", horizontal = TRUE)
                dev.off()
            }

            # plot histogram of dp and var for only identified
            a <- ggplot(y[which(!quant[, "Compound Name"] == "" & !is.na(quant[,
                "Compound Name"])), ], aes(x = y[which(!quant[, "Compound Name"] ==
                "" & !is.na(quant[, "Compound Name"])), nrow(metadata) + 1])) +
                geom_histogram() + ggtitle("Variance Distribution") + xlab("Variance") +
                ylab("Frequency") + theme(rect = element_rect(fill = "transparent"))  # histograma da variância. O ideal é que tenha pouca variância alta
            b <- ggplot(y[which(!quant[, "Compound Name"] == "" & !is.na(quant[,
                "Compound Name"])), ], aes(x = y[which(!quant[, "Compound Name"] ==
                "" & !is.na(quant[, "Compound Name"])), nrow(metadata) + 2])) +
                geom_histogram() + ggtitle("Standard Deviation") + xlab("Standard Deviation") +
                ylab("Frequency") + theme(rect = element_rect(fill = "transparent"))  # histograma do desvio padrão.

            # png
            if (".png" %in% pic_extension) {
                png("identified_variance_dp.png", units = "cm", width = 16, height = 16,
                  res = 900, bg = "NA")
                gridExtra::grid.arrange(a, b, ncol = 1)
                dev.off()
                png("identified_boxplot.png", units = "cm", width = 16, height = 16,
                  res = 900, bg = "NA")
                par(mfrow = c(2, 1))
                boxplot(y$Variance, main = "Variance", horizontal = TRUE)
                boxplot(y$sd, main = "Standard Deviation", horizontal = TRUE)
                dev.off()
            }

            # tiff
            if (".tiff" %in% pic_extension) {
                tiff("identified_variance_dp.tiff", units = "cm", width = 16,
                  height = 16, res = 900, bg = "NA")
                gridExtra::grid.arrange(a, b, ncol = 1)
                dev.off()
                tiff("identified_boxplot.tiff", units = "cm", width = 16, height = 16,
                  res = 900, bg = "NA")
                par(mfrow = c(2, 1))
                boxplot(y$Variance, main = "Variance", horizontal = TRUE)
                boxplot(y$sd, main = "Standard Deviation", horizontal = TRUE)
                dev.off()
            }

            # set up quantification data
            n <- cbind(n[, seq_len(nrow(metadata) + 5)], y[, c("Variance", "sd")])  # as everything is ok, use only first most intense ions for representative
            write.csv(n, "Normalized_quantification.csv", row.names = FALSE,
                na = "")

            res <- menu(c("Conclude normalization", "Re-do normalization", "Remove all outliers",
                "Remove spectra with variance bigger than defined"), graphics = TRUE,
                title = "Choose the best normalization")
            if (res == 1) {
                okay <- 1
            }
            if (res == 2) {
                okay <- 2
            }
            if (res == 3) {
                y_1 <- y[-which(y$Variance >= boxplot.stats(y$Variance)$stats[5]),
                  ]
                n_1 <- n[-which(y$Variance >= boxplot.stats(y$Variance)$stats[5]),
                  ]

                # png
                if (".png" %in% pic_extension) {
                  png("boxplot_withoutOutliers.png", units = "cm", width = 16,
                    height = 16, res = 900, bg = "NA")
                  par(mfrow = c(2, 1))
                  boxplot(y_1$Variance, main = "Variance", horizontal = TRUE)
                  boxplot(y_1$sd, main = "Standard Deviation", horizontal = TRUE)
                  dev.off()
                }
                # tiff
                if (".tiff" %in% pic_extension) {
                  tiff("boxplot_withoutOutliers.tiff", units = "cm", width = 16,
                    height = 16, res = 900, bg = "NA")
                  par(mfrow = c(2, 1))
                  boxplot(y_1$Variance, main = "Variance", horizontal = TRUE)
                  boxplot(y_1$sd, main = "Standard Deviation", horizontal = TRUE)
                  dev.off()
                }

                n_1 <- cbind(n_1[, seq_len(nrow(metadata) + 5)], y_1[, c("Variance",
                  "sd")])  # as everything is ok, use only first most intense ions for representative
                write.csv(n_1, "Normalized_quantification_withoutOutliers.csv",
                  row.names = FALSE, na = "")
                a <- ggplot(y_1, aes(x = y_1[, nrow(metadata) + 1])) + geom_histogram() +
                  ggtitle("Variance Distribution") + xlab("Variance") + ylab("Frequency") +
                  theme(rect = element_rect(fill = "transparent"))  # histograma da variância. O ideal é que tenha pouca variância alta
                b <- ggplot(y_1, aes(x = y_1[, nrow(metadata) + 2])) + geom_histogram() +
                  ggtitle("Standard Deviation") + xlab("Standard Deviation") +
                  ylab("Frequency") + theme(rect = element_rect(fill = "transparent"))  # histograma do desvio padrão.
                # png
                if (".png" %in% pic_extension) {
                  png("variance_dp_withoutOutliers.png", units = "cm", width = 16,
                    height = 16, res = 900, bg = "NA")
                  gridExtra::grid.arrange(a, b, ncol = 1)
                  dev.off()
                }
                # tiff
                if (".tiff" %in% pic_extension) {
                  tiff("variance_dp_withoutOutliers.tiff", units = "cm", width = 16,
                    height = 16, res = 900, bg = "NA")
                  gridExtra::grid.arrange(a, b, ncol = 1)
                  dev.off()
                }
            }

            if (res == 4) {
                th <- as.double(dlgInput("Variance threshold ", "0")$res)
                y_2 <- y[-which(y$Variance >= th), ]
                n_2 <- n[-which(y$Variance >= th), ]
                # png
                if (".png" %in% pic_extension) {
                  png(paste0("boxplot_varianceUpTo", th, ".png"), units = "cm",
                    width = 16, height = 16, res = 900, bg = "NA")
                  par(mfrow = c(2, 1))
                  boxplot(y_2$Variance, main = "Variance", horizontal = TRUE)
                  boxplot(y_2$sd, main = "Standard Deviation", horizontal = TRUE)
                  dev.off()
                }
                # tiff
                if (".tiff" %in% pic_extension) {
                  tiff(paste0("boxplot_varianceUpTo", th, ".tiff"), units = "cm",
                    width = 16, height = 16, res = 900, bg = "NA")
                  par(mfrow = c(2, 1))
                  boxplot(y_2$Variance, main = "Variance", horizontal = TRUE)
                  boxplot(y_2$sd, main = "Standard Deviation", horizontal = TRUE)
                  dev.off()
                }
                n_2 <- cbind(n_2[, seq_len(nrow(metadata) + 5)], y_2[, c("Variance",
                  "sd")])  # as everything is ok, use only first most intense ions for representative
                write.csv(n_2, paste0("Normalized_quantification_varianceUpTo",
                  th, ".csv"), row.names = FALSE, na = "")
                a <- ggplot(y_2, aes(x = y_2[, nrow(metadata) + 1])) + geom_histogram() +
                  ggtitle("Variance Distribution") + xlab("Variance") + ylab("Frequency") +
                  theme(rect = element_rect(fill = "transparent"))  # histograma da variância. O ideal é que tenha pouca variância alta
                b <- ggplot(y_2, aes(x = y_2[, nrow(metadata) + 2])) + geom_histogram() +
                  ggtitle("Standard Deviation") + xlab("Standard Deviation") +
                  ylab("Frequency") + theme(rect = element_rect(fill = "transparent"))  # histograma do desvio padrão.
                # png
                if (".png" %in% pic_extension) {
                  png("variance_dp_withoutOutliers.png", units = "cm", width = 16,
                    height = 16, res = 900, bg = "NA")
                  gridExtra::grid.arrange(a, b, ncol = 1)
                  dev.off()
                }
                # tiff
                if (".tiff" %in% pic_extension) {
                  tiff("variance_dp_withoutOutliers.tiff", units = "cm", width = 16,
                    height = 16, res = 900, bg = "NA")
                  gridExtra::grid.arrange(a, b, ncol = 1)
                  dev.off()
                }
            }
        }
        okay <- menu(c("OK, keep going", "No, re-do normalization"), graphics = TRUE,
            title = "Are the results ok?")
    }

    # set to main folder
    setwd(myDir)

    # return results
    return(n)
}
