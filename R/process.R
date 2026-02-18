#' process
#'
#' This function computes the second step of workflow (peak picking, RT correction and peak grouping) and build images for visualization of data using functions from xcms package.
#' @export
#' @param myDir Path to the directory of work.
#' @param raw_data A 'MSnExp' from MSnbase package.
#' @param metadata A matrix or data.frame with metadata information about samples. Include, at least 'sample' and 'file' columns with name of sample and its path, respectively. More information can be added in new columns, such as 'group', 'class', 'biorep' and 'tecrep'.
#' @param colors A list with colors generated from 'read_data()'.
#' @param peakMonitor Logical. Are there peak to monitor throuhout the workflow? Default to FALSE.
#' @param pictures Logical. If pictures should be plotted or not. Default to TRUE.
#' @param example Logical. If is example, pop-ups won't appear and images won't be generated.
#' @param ions List with sublist mz = mz (numeric) of the monitored ion and rt = retention time of monitored ion (numeric). To the 'rt' will be added and subtracted 5 seconds. Default to null.
#' @param filter Numeric. Intensity threshold for the peak detection. Default to NULL. When NULL, the user will be asked for a number. Set filter = 0 for no intensity filtering.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param group Character. Name from 'metadata' column names to group the samples. Default to NULL.
#' @return A 'xcmsSet' object with detected, grouped and filled peaks with retention time corrected.
#' @importFrom grDevices dev.off pdf png tiff
#' @importFrom graphics boxplot grid legend par text
#' @importFrom methods as new
#' @importFrom stringr str_c
#' @importFrom svDialogs dlgInput dlg_list
#' @importFrom xcms MatchedFilterParam findChromPeaks refineChromPeaks FilterIntensityParam adjustRtime ObiwarpParam groupChromPeaks PeakDensityParam fillChromPeaks ChromPeakAreaParam plotChromPeakImage chromPeaks plotAdjustedRtime chromatogram fillPeaks
#' @importFrom utils choose.dir memory.limit menu read.csv select.list write.csv write.table
#' @examples
#' \donttest{
#' \dontrun{
#' load(system.file("extdata", "raw_data.RData", package = "PipMet"))
#' load(system.file("extdata", "metadata.RData", package = "PipMet"))
#' load(system.file("extdata", "colors.RData", package = "PipMet"))
#' xdata4 <- process(raw_data, metadata, myDir = NULL, colors, example = TRUE, group = "group")
#' }
#' }
process <- function(raw_data, metadata, myDir, colors, peakMonitor = FALSE, ions = NULL,
                    pictures = TRUE, example = FALSE, filter = NULL, pic_extension = c(
                      ".tiff",
                      ".png"
                    ), group = NULL) {
  message(str_c("Processing ", length(raw_data$files), " files..."))
  # peak picking
  message("Detecting peaks...")
  mfp <- MatchedFilterParam(
    fwhm = 5, binSize = 0.5, steps = 2, mzdiff = 0.5,
    snthresh = 2, max = 500
  )
  xdata <- findChromPeaks(raw_data, param = mfp)
  if (!example == TRUE) {
    save(xdata, file = "xdata.RData")
  }

  if (example == FALSE) {
    if (!is.null(filter)) {
      message("Filtering peaks...")
      xdata <- refineChromPeaks(xdata, param = FilterIntensityParam(
        threshold = as.integer(filter),
        nValues = 1, value = "maxo"
      ))
    } else {
      filt <- menu(c("Yes", "No"), graphics = TRUE, title = "Apply intensity filter?")
      if (filt == 1) {
        filter <- dlgInput("Intensity threshold ", "0")$res
        message("Filtering peaks...")
        xdata <- refineChromPeaks(xdata, param = FilterIntensityParam(
          threshold = as.integer(filter),
          nValues = 1, value = "maxo"
        ))
      }
    }
  }

  # check if there is more than one file
  if (nrow(xdata@phenoData@data) > 1) {
    # retention time correction
    message("Adjusting retention time along samples...")
    xdata2 <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
    if (!example == TRUE) {
      save(xdata2, file = "xdata2.RData")
    }

    if (example == FALSE) {
      if (is.null(group)) {
        # ask condition to compare from the metadata table
        group <- dlg_list(colnames(metadata), multiple = FALSE, title = "Condition to group from:")$res
      }
      # group <- dlg_list(colnames(metadata), multiple = FALSE, title
      # = 'Choose conditions (from metadata table) to group samples:
      # ')$res
    } else {
      group <- "group"
    }

    # grouping peaks
    message("Grouping peaks along samples...")
    xdata3 <- groupChromPeaks(xdata2, param = PeakDensityParam(sampleGroups = metadata[
      ,
      group
    ], bw = 0.5, minSamples = 1, maxFeatures = 500, minFraction = 0.4))
    if (!example == TRUE) {
      save(xdata3, file = "xdata3.RData")
    }

    # imagens dos dados em processamento e pós processamento (padrões,
    # cromatogramas de íon extraído)
    if (exists("xdata3") & pictures == TRUE & example == FALSE) {
      # create and set folder for images
      if (!dir.exists("peakProcessing_results")) {
        dir.create("peakProcessing_results")
      }
      setwd("peakProcessing_results")

      # heatmap of identified peaks per region of chromatogram tiff
      if (".tiff" %in% pic_extension) {
        tiff("plotChromPeakImage.tiff",
          units = "cm", width = 16, height = 16,
          res = 900, bg = "NA"
        )
        par(mar = c(5, 9, 4, 1) + 0.1)
        plotChromPeakImage(xdata3)
        dev.off()
      }
      # png
      if (".png" %in% pic_extension) {
        png("plotChromPeakImage.png",
          units = "cm", width = 16, height = 16,
          res = 900, bg = "NA"
        )
        par(mar = c(5, 9, 4, 1) + 0.1)
        plotChromPeakImage(xdata3)
        dev.off()
      }

      # boxplot of log2 intensities per sample
      ints <- split(log2(chromPeaks(xdata3)[, "into"]), f = chromPeaks(xdata3)[,"sample"])
      names(ints) <- metadata$sample
      for (i in seq_len(length(colors))) {
        # tiff
        if (".tiff" %in% pic_extension) {
          tiff(paste0(names(colors)[i], "_boxplotLog2Postprocessed.tiff"),
            units = "cm", width = 16, height = 16, res = 900, bg = "NA"
          )
          par(mar = c(7, 5, 3, 1) + 0.1, cex.axis = 1, xpd = TRUE)
          boxplot(ints,
            varwidth = TRUE, col = colors[[i]][[2]][colors[[i]][[1]]],
            ylab = expression(log[2] ~ intensity), main = "Peak intensities",
            las = 3, xaxt = "n"
          )
          grid(nx = NA, ny = NULL)
          text(seq_along(metadata$sample), par("usr")[3],
            labels = metadata$sample,
            srt = 45, adj = c(1.1, 1.1), xpd = TRUE, cex = 0.7
          )
          legend("right",
            inset = c(-0.23, 0), legend = names(colors[[i]][[2]]), ,
            box.lty = 0, col = colors[[i]][[2]], fill = colors[[i]][[2]],
            bg = "transparent"
          )

          dev.off()
        }
        # png
        if (".png" %in% pic_extension) {
          png(paste0(names(colors)[i], "_boxplotLog2Postprocessed.png"),
            units = "cm", width = 16, height = 16, res = 900, bg = "NA"
          )
          par(mar = c(7, 5, 3, 8) + 0.1, cex.axis = 1, xpd = TRUE)
          boxplot(ints,
            varwidth = TRUE, col = colors[[i]][[2]][colors[[i]][[1]]],
            ylab = expression(log[2] ~ intensity), main = "Peak intensities",
            las = 3, xaxt = "n"
          )
          grid(nx = NA, ny = NULL)
          text(seq_along(metadata$sample), par("usr")[3],
            labels = metadata$sample,
            srt = 45, adj = c(1.1, 1.1), xpd = TRUE, cex = 0.65
          )
          legend("right",
            inset = c(-0.23, 0), legend = names(colors[[i]][[2]]), ,
            box.lty = 0, col = colors[[i]][[2]], fill = colors[[i]][[2]],
            bg = "transparent"
          )
          dev.off()
        }
      }

      # chromatogram postprocessed
      bpc_after <- chromatogram(xdata3, aggregationFun = "max", include = "none")
      for (i in seq_len(length(colors))) {
        # tiff
        if (".tiff" %in% pic_extension) {
          tiff(paste0(names(colors)[i], "_postprocessedChromatogram.tiff"),
            units = "cm", width = 16, height = 16, res = 900, bg = "NA"
          )
          par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
          plot(bpc_after, col = colors[[i]][[2]][colors[[i]][[1]]])
          legend("topright",
            legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
            fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent"
          )
          plotAdjustedRtime(xdata3, col = colors[[i]][[2]][colors[[i]][[1]]])
          legend("bottomright",
            legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
            fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent"
          )
          dev.off()
        }
        # png
        if (".png" %in% pic_extension) {
          png(paste0(names(colors)[i], "_postprocessedChromatogram.png"),
            units = "cm", width = 16, height = 16, res = 900, bg = "NA"
          )
          par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
          plot(bpc_after, col = colors[[i]][[2]][colors[[i]][[1]]])
          legend("topright",
            legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
            fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent"
          )
          plotAdjustedRtime(xdata3, col = colors[[i]][[2]][colors[[i]][[1]]])
          legend("bottomright",
            legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
            fill = colors[[i]][[2]], box.lty = 0, cex = 1, bg = "transparent"
          )
          dev.off()
        }
      }

      # cromatograma de íons extraído
      if (peakMonitor == TRUE | peakMonitor == "Yes") {
        if (!dir.exists("Monitoring ions") == TRUE) {
          dir.create("Monitoring ions")
        }
        setwd("Monitoring ions")

        for (ii in seq_len(length(ions))) {
          crom <- chromatogram(xdata3, rt = c(
            ions[[ii]][["rt"]] - 5,
            ions[[ii]][["rt"]] + 5
          ), mz = ions[[ii]][["mz"]], include = "none")
          for (i in seq_len(length(colors))) {
            # tiff
            if (".tiff" %in% pic_extension) {
              tiff(paste0(names(colors)[i], "_", ii, "_postPross_EIC.tiff"),
                units = "cm", width = 16, height = 16, res = 900, bg = "NA"
              )
              plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
              legend("right",
                legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent"
              )
              dev.off()
            }
            # png
            if (".png" %in% pic_extension) {
              png(paste0(names(colors)[i], "_", ii, "_postPross_EIC.png"),
                units = "cm", width = 16, height = 16, res = 900, bg = "NA"
              )
              plot(crom, col = colors[[i]][[2]][colors[[i]][[1]]])
              legend("right",
                legend = names(colors[[i]][[2]]), col = colors[[i]][[2]],
                fill = colors[[i]][[2]], box.lty = 0, cex = 0.8, bg = "transparent"
              )
              dev.off()
            }
          }
        }
      }

      # set to main folder
      setwd(myDir)
    }
  }

  if (exists("xdata3")) {
    xdata4 <- as(xdata3, "xcmsSet")
    if ("group" %in% colnames(metadata)) {
      sampclass(xdata4) <- metadata$group
    } else {
      if ("class" %in% colnames(metadata)) {
        sampclass(xdata4) <- metadata$class
      } else {
        sampclass(xdata4) <- NA
      }
    }
    message("Filling missing peaks...")
    xdata4 <- fillPeaks(xdata4)
    if (!example == TRUE) {
      save(xdata4, file = "xdata4.RData")
    }
  } else {
    xdata4 <- as(xdata, "xcmsSet")
  }

  # return results
  return(list(xdata4 = xdata4, group = group))
}
