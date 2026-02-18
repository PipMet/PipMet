#' read_data
#'
#' This function creates a sheet for the user to fill with the experiment design in order to create categories for the files to be processed. Also creates images of preprocessed data.
#' @keywords metadata
#' @export
#' @return A list of list of colors, metadata table, a OnDiskMSnExp object and the path to the directory of work.
#' @param ions List with sublist mz = mz (numeric) of the monitored ion and rt = retention time of monitored ion (numeric). To the 'rt' will be added and subtracted 5 seconds. Default to NULL.
#' @param myDir Path to working directory
#' @param sample_dir Path to sample directory.
#' @param metadata Path to .csv file or an R data.frame object containing metadata. At least 'sample' and 'file' columns must be included.
#' @param extension Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'.
#' @param pictures Logical. If pictures should be plotted or not. Default to TRUE.
#' @param example Logical. If example = TRUE, the metadata and other needed files will be loaded from package files.
#' @param peakMonitor Logical. Are there peak to monitor throuhout the workflow? Default to FALSE.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param sample_names Character. Name of metadata column to name the samples. Default to NULL. If NULL, the user will be asked.
#' @importFrom grDevices dev.off png tiff rainbow
#' @importFrom graphics boxplot legend par
#' @importFrom methods as new
#' @importFrom stringr str_c
#' @importFrom utils choose.dir menu read.csv write.csv write.table choose.files
#' @importFrom svDialogs dlgInput dlg_message dlg_list
#' @importFrom fritools is_path
#' @import xcms
#' @import MSnbase
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \donttest{
#' \dontrun{
#' read <- read_data(
#'   pictures = FALSE,
#'   example = TRUE
#' )
#' colors <- read[[1]]
#' metadata <- read[[2]]
#' raw_data <- read[[3]]
#' myDir <- read[[4]]
#' rm(read)
#' 
#' }}
read_data <- function(peakMonitor = NULL, ions = NULL, myDir = NULL, sample_dir = NULL,
    metadata = NULL, extension = NULL, pictures = TRUE, example = FALSE, pic_extension = c(".tiff",
        ".png"), sample_names = NULL) {

    if (!example == TRUE) {

        # ask user about samples and folders path and create a new folder
        # named after a 'Project'
        if (is.null(sample_dir) | missing(sample_dir)) {
            sample_dir <- choose.dir(default = getwd(), caption = "Please, select the Samples directory, should be C:/Users/_/Samples")
        }
        setwd(sample_dir)
        sample_dir <- getwd()
        setwd("~/")

        if (is.null(myDir) | missing(myDir)) {
            myDir <- dlgInput("Name your project", Sys.info()["user"])$res
            dir.create(myDir, showWarnings = FALSE)
        }
        setwd(myDir)
        myDir <- getwd()

        # file extension
        if (is.null(metadata)) {
            if (is.null(extension)) {
                extension <- dlg_list(c(".mzML", ".mzXML"), multiple = TRUE,
                  title = "Files extension:")$res
            }
        }

        # set up metadata table
        if (!is.null(metadata)) {
            metadata <- read.csv(metadata, na.string = c("NA", ""), colClasses = "character",
                sep = ",")
        }

        if (is.null(metadata)) {
            metd <- dlg_message("Metadata table already exists?", "yesno")$res
            if (metd == "yes") {
                metadata <- read.csv(choose.files(), na.string = c("NA", ""),
                  colClasses = "character", sep = ",")
                if (!"file" %in% colnames(metadata)) {
                  metadata$file <- paste0(sample_dir, "/", metadata$sample, extension)
                  setwd(myDir)
                } else {
                  for (i in seq_len(nrow(metadata))) {
                    # check if they are just filenames or filepath (they
                    # need to be path)
                    if (!is_path(metadata$file[i])) {
                      setwd(sample_dir)
                      metadata$file[i] <- file_path_as_absolute(basename(metadata$file[i]))
                    } else {
                      metadata$file[i] <- paste0(sample_dir, "/", basename(metadata$file[i]))
                    }
                  }
                  setwd(myDir)
                }
            } else {
                files <- list.files(sample_dir, full.names = FALSE, pattern = extension,
                  recursive = TRUE)
                metadata <- matrix(nrow = length(files), ncol = 6)
                colnames(metadata) <- c("sample", "group", "class", "tec_rep",
                  "bio_rep", "file")
                metadata[, "file"] <- files
                metadata[, "sample"] <- sub(basename(files), pattern = extension,
                  replacement = "", fixed = TRUE)
                setwd(sample_dir)
                for (i in seq_len(nrow(metadata))) {
                  metadata[i, "file"] <- file_path_as_absolute(metadata[i, "file"])
                }
                setwd(myDir)
                write.csv(metadata, "metadata.csv", row.names = FALSE)
                dlg_message("A file 'metadata.csv' was created in you directory. Fill the sheet before continuing. In 'class' column, describe your samples as 'Sample', 'QC', 'Blank' or 'Pool'. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.",
                  type = "ok")
                while (file.exists("metadata.csv") == FALSE) {
                  dlg_message("A file 'metadata.csv' was created in you directory. Fill the sheet before continuing. In 'class' column, describe your samples as 'Sample', 'QC', 'Blank' or 'Pool'. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.",
                    type = "ok")
                  metadata <- read.csv("metadata.csv", na.string = c("NA", ""),
                    colClasses = "character", sep = ",")
                }
                metadata <- read.csv("metadata.csv", na.string = c("NA", ""),
                  colClasses = "character", sep = ",")
            }
            if (!sum((is.na(metadata))) == 0) {
                dlg_message("It seems to existx empty columns/rows in you metadata file. Please, delete and press 'OK'.")$res
                metadata <- read.csv(choose.files(), na.string = c("NA", ""),
                  colClasses = "character", sep = ",")
                if (length(colnames(metadata))<=1) {
                    metadata <- read.csv(choose.files(), na.string = c("NA", ""),
                        colClasses = "character", sep = ";")}
            }
        }
    } else {
        # if example == TRUE
        metadata <- read.csv(system.file("extdata", "metadata.csv", package = "PipMet"),
            na.string = c("NA", ""), colClasses = "character", sep = ",", dec = ".")
        extension <- ".mzXML"
        for (i in seq_len(nrow(metadata))) {
            metadata$file[i] <- system.file("extdata", metadata$file[i], package = "PipMet")
        }
        sample_dir <- system.file("extdata", package = "PipMet")
        myDir <- dlgInput("Name your project", "PipMet_example")$res
        if (!dir.exists(myDir) == TRUE) {
            dir.create(myDir, showWarnings = FALSE)
        }
        setwd(myDir)
        myDir <- getwd()
    }

    # if there is only a single file
    if (!nrow(metadata) == 1) {
        # create 'metadata$all' 1 and 2 for identification
        x <- colnames(metadata)
        for (i in c("file", "all", "all2")) {
            if (i %in% x) {
                x <- x[-which(x == i)]
            }
        }

        for (i in x) {
            if (length(unique(metadata[, i])) > 10) {
                x <- x[-which(x == i)]
            }
        }

        # colors for each column in metadata except 'sample' and 'tec_rep'
        colors <- vector(mode = "list", length = length(x))
        names(colors) <- x
        for (i in seq_len(length(x))) {
            if (length(unique(metadata[, x[[i]]])) <= 9) {
                colors[[i]] <- list(metadata[, x[i]], paste0(RColorBrewer::brewer.pal(length(unique(metadata[,
                  x[[i]]])), "Set1")[seq_along(unique(metadata[, x[[i]]]))], "60"))
            } else {
                colors[[i]] <- list(metadata[, x[i]], rainbow(length(unique(metadata[,
                  x[[i]]]))))
            }
            names(colors[[i]][[2]]) <- c(unique(metadata[, x[i]]))
            names(colors[[i]]) <- c(x[[i]], paste0(x[[i]], "_colors"))
        }
    }

    # read data into R
    message(str_c("Reading ", nrow(metadata), " files..."))
    raw_data <- readMSData(metadata$file, pdata = new("AnnotatedDataFrame",
        metadata), mode = "onDisk")

    if (pictures == TRUE) {
        if (is.null(sample_names)) {
            sample_names <- dlg_list(names(metadata), multiple = FALSE, title = "Name sample as:")$res
        }
        # images of pre-processing
        if (!dir.exists("Visualization_results") == TRUE) {
            dir.create("Visualization_results")
        }
        setwd("Visualization_results")

        # get info to plot
        bpc <- chromatogram(raw_data, aggregationFun = "max")
        tic <- chromatogram(raw_data, aggregationFun = "sum")

        # chromatograms
        for (i in seq_len(length(colors))) {
            # tiff
            if (".tiff" %in% pic_extension) {
                tiff(paste0(names(colors)[i], "_chromatograms.tiff"), units = "cm",
                  width = 16, height = 16, res = 900, bg = "NA")
                par(mfrow = c(2, 1))
                plot(bpc, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Base Peak Chromatogram")
                legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                  fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
                plot(tic, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Total Ion Current Chromatogram")
                legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                  fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
                dev.off()
            }
            # png
            if (".png" %in% pic_extension) {
                png(paste0(names(colors)[i], "_chromatograms.png"), units = "cm",
                  width = 16, height = 16, res = 900, bg = "NA")
                par(mfrow = c(2, 1))
                plot(bpc, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Base Peak Chromatogram")
                legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                  fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
                plot(tic, col = colors[[i]][[2]][colors[[i]][[1]]], main = "Total Ion Current Chromatogram")
                legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                  fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
                dev.off()
            }
        }

        # boxplot of total ion current - FIX tic_por_arquivo <-
        # split(tic(raw_data), f = fromFile(raw_data)) for (i in
        # seq_len(length(colors))) { tic_por_arquivo <-
        # split(tic(raw_data), f = fromFile(raw_data)) # tiff
        # tiff(paste0(names(colors)[i], '_ticBoxplot.tiff'), units = 'cm',
        # width = 16, height = 16, res = 900, bg = 'NA')
        # boxplot(tic_por_arquivo, col =
        # colors[[i]][[2]][colors[[i]][[1]]], ylab = 'intensity', xlab =
        # 'sample', main = 'Total ion current') dev.off() # png
        # png(paste0(names(colors)[i], '_ticBoxplot.png'), units = 'cm',
        # width = 16, height = 16, res = 900, bg = 'NA')
        # boxplot(tic_por_arquivo, col =
        # colors[[i]][[2]][colors[[i]][[1]]], ylab = 'intensity', xlab =
        # 'sample', main = 'Total ion current') dev.off() }
        # rm(tic_por_arquivo)

        if (!nrow(metadata) == 1) {
            # cluster
            tic_bin <- bin(tic, binSize = 1)
            cl <- do.call(cbind, lapply(tic_bin, intensity))
            cl[cl == 0] <- NA
            cormat <- cor(log2(cl), use = "pairwise.complete.obs")
            colnames(cormat) <- rownames(cormat) <- metadata[, sample_names]
            # colnames(cormat) <- rownames(cormat) <- metadata[, 'sample']
            # for each set of colors (conditions of experiment)
            for (i in seq_len(length(colors))) {
                ann <- data.frame(colors[[i]][[1]])
                colnames(ann) <- names(colors)[i]
                rownames(ann) <- metadata[, "sample"]
                ant <- list(colors[[i]][[2]])
                names(ant) <- names(colors)[i]
                # tiff
                if (".tiff" %in% pic_extension) {
                  tiff(paste0(names(colors)[i], "_", "sample", "_", "_cluster.tiff"),
                    units = "cm", width = 16, height = 16, res = 900, bg = "NA")
                  pheatmap(cormat, annotation = ann, annotation_colors = ant,
                    border_color = "NA", cluster_rows = FALSE, )
                  dev.off()
                }
                # png
                if (".png" %in% pic_extension) {
                  png(paste0(names(colors)[i], "_", "sample", "_", "_cluster.png"),
                    units = "cm", width = 16, height = 16, res = 900, bg = "NA")
                  pheatmap(cormat, annotation = ann, annotation_colors = ant,
                    border_color = "NA", cluster_rows = FALSE, )
                  dev.off()
                }
            }
            rm(tic_bin)
        }

        # extracted ion chromatogram based on mz and rt asked previously by
        # user
        if (peakMonitor == TRUE | peakMonitor == "Yes") {
            if (!dir.exists("Monitoring ions") == TRUE) {
                dir.create("Monitoring ions")
            }
            setwd("Monitoring ions")

            for (ii in seq_len(length(ions))) {
                crom <- chromatogram(raw_data, rt = c(as.numeric(ions[[ii]][["rt"]] -
                  5), as.numeric(ions[[ii]][["rt"]] + 5)), mz = as.numeric(ions[[ii]][["mz"]]))
                for (i in seq_len(length(colors))) {
                  # tiff
                  if (".tiff" %in% pic_extension) {
                    tiff(paste0(names(colors)[i], "_", ii, "_prePross_EIC.tiff"),
                      units = "cm", width = 16, height = 16, res = 900, bg = "NA")
                    plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
                    legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                      fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
                    dev.off()
                  }

                  # png
                  if (".png" %in% pic_extension) {
                    png(paste0(names(colors)[i], "_", ii, "_prePross_EIC.png"),
                      units = "cm", width = 16, height = 16, res = 900, bg = "NA")
                    plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
                    legend("right", legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                      fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent")
                    dev.off()
                  }
                }
            }
        }
    }


    # return to main folder
    setwd(myDir)

    # return results
    return(list(colors = colors, metadata = metadata, raw_data = raw_data, myDir = myDir))
}
