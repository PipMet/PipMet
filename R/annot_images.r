#' Annotation images
#'
#' This function creates network, hierarchy and heatmap plots of the spectra (annotated or not).
#' @export
#' @param example Default to FALSE.
#' @param pslist List of spectra.
#' @param myDir Path to the directory of work.
#' @param pictures Logical. If pictures should be plotted or not. Default to TRUE.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param pre_anno Path to pre_anno.csv annotation file or pre_anno table. If NULL, the pre_anno.csv in the folder will be read. Default to NULL.
#' @return A list with (1) 'pseudoespectrum' object annotated resulted from 'addAnnotations' (CluMSID package) and (2) a table with annotation.
#' @importFrom grDevices dev.off pdf png tiff
#' @importFrom graphics par
#' @importFrom utils read.csv
#' @importFrom fritools is_path
#' @importFrom CluMSID addAnnotations distanceMatrix networkplot HCplot
#' @examples
#' \donttest{
#' \dontrun{
#' load(system.file("extdata", "pslist.RData", package = "PipMet"))
#' annot <- annot_images(pslist,
#'   example = TRUE,
#'   pictures = FALSE,
#'   pre_anno = system.file("extdata", "pre_anno.csv",
#'     package = "PipMet"
#'   )
#' )
#' }
#' }
annot_images <- function(pslist = NULL, myDir = "~/", pictures = TRUE, pic_extension = c(
                           ".tiff",
                           ".png"
                         ), pre_anno = NULL, example = FALSE) {
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
      r <- read.csv("pre_anno.csv",
        sep = ";", na.string = c("NA", ""),
        dec = ","
      )
    }
  }
  if (sum(is.na(r$id)) > 0) {
    r <- r[-which(is.na(r$id)), ]
  }
  r[is.na(r)] <- ""

  ### add annotations from 'pre_anno.csv' file (if there is annotation)

  if (sum(is.na(r$annotation)) == nrow(r)) {
    apslist <- pslist
  }
  apslist <- addAnnotations(featlist = pslist, annolist = r)

  if (pictures == TRUE & nrow(r) > 5 & !example == TRUE) {
    pseudodistmat <- distanceMatrix(apslist, mz_tolerance = 0.02) ### calculates distance matrix; takes a while

    # creates a folder for images
    if (!dir.exists("Statistics") == TRUE) {
      dir.create("Statistics")
    }
    setwd("Statistics")

    # network plot tiff tiff('network_plot_0.7.tiff', units = 'cm',
    # width = 16, height = 16, res = 900, bg = 'NA')
    # networkplot(pseudodistmat, highlight_annotated = TRUE,
    # show_labels = TRUE, exclude_singletons = TRUE, min_similarity =
    # 0.7) dev.off() png png('network_plot_0.7.png', units = 'cm',
    # width = 16, height = 16, res = 900, bg = 'NA')
    # networkplot(pseudodistmat, highlight_annotated = TRUE,
    # show_labels = TRUE, exclude_singletons = TRUE, min_similarity =
    # 0.7) dev.off()

    # tiff
    if (".tiff" %in% pic_extension) {
      # hierarchy plot
      tiff("hierarchy_plot.tiff",
        units = "cm", width = 16, height = 16,
        res = 900, bg = "NA"
      )
      par(mar = c(4.5, 4.2, 1, 0.5))
      HCplot(pseudodistmat, h = 0.7, cex = 0.2)
      dev.off()
      # heatmap plot
      tiff("heatmap.tiff",
        units = "cm", width = 16, height = 16, res = 900,
        bg = "NA"
      )
      par(mar = c(4.5, 4.2, 1, 0.5))
      HCplot(pseudodistmat,
        type = "heatmap", cexRow = 0.2, cexCol = 0.24,
        margins = c(7, 7)
      )
      dev.off()
    }

    # png
    if (".png" %in% pic_extension) {
      # hierarchy plot
      png("hierarchy_plot.png",
        units = "cm", width = 16, height = 16,
        res = 900, bg = "NA"
      )
      par(mar = c(4.5, 4.2, 1, 0.5))
      HCplot(pseudodistmat, h = 0.7, cex = 0.2)
      dev.off()
      # heatmap plot
      png("heatmap.png",
        units = "cm", width = 16, height = 16, res = 900,
        bg = "NA"
      )
      par(mar = c(4.5, 4.2, 1, 0.5))
      HCplot(pseudodistmat,
        type = "heatmap", cexRow = 0.2, cexCol = 0.2,
        margins = c(7, 7)
      )
      dev.off()
    }
  }

  # set to main folder
  if (!example == TRUE) {
    setwd(myDir)
  }

  # return results
  return(list(apslist = apslist, r = r))
}
