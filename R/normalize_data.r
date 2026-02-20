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
#' Normalize data
#'

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
            okay_merge <- "yes"
            while (okay_merge == "yes") {
                z_merge <- dlg_list(sort(quant[which(!is.na(quant[, "Compound Name"]) &
                  !quant[, "Compound Name"] == ""), "Compound Name"]), multiple = TRUE,
                  title = "Merge")$res
                if (length(z_merge) >= 2) {
                  y_merge <- vector()
                  for (i in seq_along(unique(z_merge))) {
                    y_merge <- rbind(y_merge, quant[which(quant[, "Compound Name"] == unique(z_merge)[i]),
                      ])
                  }
                  x_merge <- rbind(y_merge[which(y_merge[, "id"] == sort(as.integer(y_merge[, "id"]))[1]),
                    ])
                  for (i in 6:ncol(quant)) {
                    x_merge[1, i] <- sum(unlist(lapply(y_merge[, i], as.numeric)))
                  }
                  quant[which(quant[, "id"] == sort(as.integer(y_merge[, "id"]))[1]),
                    ] <- x_merge[1, ]
                  for (i in 2:length(z_merge)) {
                    quant <- quant[-which(quant[, "id"] == sort(as.integer(y_merge[,
                      "id"]))[i]), ]
                  }
                  okay_merge <- dlg_message("Repeat?", type = "yesno")$res
                } else {
                  okay_merge <- "no"
                }
            }
        }

        # remove any compound identified?
        if (is.null(removeCompounds)) {
            a_rem <- dlg_message("Remove compounds?", type = "yesno")$res
        } else {
            if (removeCompounds == TRUE) {
                a_rem <- "yes"
            }
            if (removeCompounds == FALSE) {
                a_rem <- "no"
            }
        }
        if (a_rem == "yes") {
            okay_rem <- "yes"
            while (okay_rem == "yes") {
                z_rem <- dlg_list(sort(quant[which(!is.na(quant[, "Compound Name"]) &
                  !quant[, "Compound Name"] == ""), "Compound Name"]), multiple = TRUE,
                  title = "Remove")$res
                if (length(z_rem) > 0) {
                    quant <- quant[-which(quant[, "Compound Name"] %in% z_rem), ]
                }
                okay_rem <- dlg_message("Repeat?", type = "yesno")$res
            }
        }

        # ----------------------------------------------------------------------
        # BLOCO INVERTIDO: NORMALIZAÇÃO POR OUTRO CRITÉRIO (MASSA, ETC) AGORA VEM ANTES
        # ----------------------------------------------------------------------
        write.csv(quant, "NotNormalized_quantification.csv", row.names = FALSE, na = "")
        g_mass <- dlg_message("Normalize by mass/numer of cells/etc?", "yesno")$res
        if (g_mass == "yes") {
            dlg_message("Normalize the \"NotNormalized_quantification.csv\" file as you wish and press \"ok\" to continue.", "ok")$res
            quant <- read.csv("NotNormalized_quantification.csv", na = "", check.names = FALSE)
        }
        # ----------------------------------------------------------------------

        if (nrow(quant) > 1 & nrow(metadata) > 1) {
            # get the quantification data
            mat <- quant[, 6:ncol(quant)]
            mat <- apply(mat, c(1, 2), FUN = as.numeric)

            # Normalyzer: evaluation and picking of normalization method
            x_norm <- as.data.frame(mat)
            x_norm[x_norm == 0] <- NA
            x_norm[x_norm > 0 & x_norm < 1] <- 0
            
            # Use only the first ion columns for NormalyzerDE evaluation
            write.table(x_norm[, seq_len(nrow(metadata))], "data.tsv", sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
            write.table(metadata, "design.tsv", sep = "\t", col.names = TRUE,
                row.names = FALSE, quote = FALSE)

            designFp <- file_path_as_absolute("design.tsv")
            dataFp <- file_path_as_absolute("data.tsv")

            if (is.null(group)) {
                group <- dlg_list(colnames(metadata), multiple = FALSE, title = "Condition to group from:")$res
            }
            
            if (nrow(mat) > 50) {
                normalyzer(jobName = "Normalyzer_results", designPath = designFp,
                  dataPath = dataFp, outputDir = getwd(), sampleColName = "sample",
                  groupColName = group, requireReplicates = FALSE)
            } else {
                try(normalyzer(jobName = "Normalyzer_results", designPath = designFp,
                  dataPath = dataFp, outputDir = getwd(), sampleColName = "sample",
                  groupColName = group, requireReplicates = FALSE, sampleAbundThres = nrow(mat)))
                try(dev.off(), silent = TRUE)
            }

            fill <- list.files(paste0(getwd(), "/Normalyzer_results"), full.names = TRUE,
                pattern = "-normalized.txt", recursive = TRUE)
            norms <- sub("-normalized.txt", replacement = "", fixed = TRUE, x = basename(fill))

            bestNormMat <- menu(norms, graphics = TRUE, title = "Choose the best normalization")

            mat[mat == 0] <- NA
            mat[mat > 0 & mat < 1] <- 0
            
            # Apply chosen method
            method <- norms[bestNormMat]
            if (method == "CycLoess") mat <- performCyclicLoessNormalization(mat)
            if (method == "GI") mat <- globalIntensityNormalization(mat)
            if (method == "log2") mat <- log2(mat)
            if (method == "mean") mat <- meanNormalization(mat)
            if (method == "median") mat <- medianNormalization(mat)
            if (method == "Quantile") mat <- performQuantileNormalization(mat)
            if (method == "RLR") mat <- performGlobalRLRNormalization(mat)
            if (method == "VSN") mat <- performVSNNormalization(mat)
            
            n <- as.data.frame(cbind(quant[, seq_len(5)], mat))

            # Statistics (Variance and SD)
            y_stats <- data.frame(matrix(nrow = nrow(n), ncol = nrow(metadata)))
            for (ii in 6:(5 + nrow(metadata))) {
                for (i in seq_len(nrow(n))) {
                  y_stats[i, (ii - 5)] <- as.double(n[i, ii]) / as.double(n[i, (ii + nrow(metadata))])
                }
            }
            y_stats[, nrow(metadata) + 1] <- apply(y_stats[, seq_len(nrow(metadata))], MARGIN = 1, FUN = var, na.rm = TRUE)
            y_stats[, nrow(metadata) + 2] <- apply(y_stats[, seq_len(nrow(metadata))], MARGIN = 1, FUN = sd, na.rm = TRUE)
            rownames(y_stats) <- n$id
            colnames(y_stats) <- c(metadata$sample, "Variance", "sd")

            # Plots
            a_plot <- ggplot(y_stats, aes(x = y_stats[, nrow(metadata) + 1])) + geom_histogram() +
                ggtitle("Variance Distribution") + xlab("Variance") + ylab("Frequency") +
                theme(rect = element_rect(fill = "transparent"))
            b_plot <- ggplot(y_stats, aes(x = y_stats[, nrow(metadata) + 2])) + geom_histogram() +
                ggtitle("Standard Deviation") + xlab("Standard Deviation") +
                ylab("Frequency") + theme(rect = element_rect(fill = "transparent"))

            if (".png" %in% pic_extension) {
                png("variance_dp.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
                gridExtra::grid.arrange(a_plot, b_plot, ncol = 1); dev.off()
                png("boxplot.png", units = "cm", width = 16, height = 16, res = 900, bg = "NA")
                par(mfrow = c(2, 1))
                boxplot(y_stats$Variance, main = "Variance", horizontal = TRUE)
                boxplot(y_stats$sd, main = "Standard Deviation", horizontal = TRUE); dev.off()
            }

            # Final results setup
            n <- cbind(n[, seq_len(nrow(metadata) + 5)], y_stats[, c("Variance", "sd")])
            write.csv(n, "Normalized_quantification.csv", row.names = FALSE, na = "")

            res <- menu(c("Conclude normalization", "Re-do normalization", "Remove all outliers",
                "Remove spectra with variance bigger than defined"), graphics = TRUE,
                title = "Choose the best normalization")
            
            if (res == 1) okay <- 1
            if (res == 2) okay <- 2
            if (res == 3) {
                # Logic for outliers...
                okay <- 1
            }
            if (res == 4) {
                # Logic for threshold...
                okay <- 1
            }
        }
        okay <- menu(c("OK, keep going", "No, re-do normalization"), graphics = TRUE,
            title = "Are the results ok?")
    }

    setwd(myDir)
    return(n)
}