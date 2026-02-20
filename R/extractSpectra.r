#' Extract Spectra
#'
#' This is a modified function from CluMSID package to extract the spectra from a CAMERA object output.
#' Not meant to be called by the user.
#' @return None.
#' @param x A 'CAMERA' xsAnnotate object.
#' @param min_peaks Numeric. Minimal number of peaks a spectrum must have to be considered a viable spectrum. Default to 5.
#' @importFrom methods as new is
#' @importFrom Biobase rowMax
#' @importFrom stats median
#' @importClassesFrom CluMSID pseudospectrum
#' @keywords pseudoespectra CAMERA spectra spectrum.
#' @examples
#' \donttest{
#' \dontrun{
#' load(system.file("extdata", "anIC.RData", package = "PipMet"))
#' spectra <- extractSpectra(anIC)
#' }
#' }
#'
extractSpectra <- function(x, min_peaks = 5) {
  # intensidade: integral do pico
  if (isTRUE(is(x, "xsAnnotate")) == TRUE) {
    if (isTRUE(sum(grepl("X", colnames(x@groupInfo)) == TRUE) == 0) == TRUE) {
      # for database, theres is no sample column in x@groupInfo
      sv <- which(colnames(x@groupInfo) == "into")
    } else {
      # search for columns in x@groupInfo correspondent to each
      # sample
      temp <- list()
      for (j in seq_along(rownames(x@xcmsSet@phenoData))) {
        temp[[j]] <- grep(paste0(
          "^X", rownames(x@xcmsSet@phenoData)[j],
          "$"
        ), colnames(x@groupInfo))
      }
      sv <- unique(unlist(temp))
    }

    # remove NAs
    x@groupInfo[, sv][is.na(x@groupInfo[, sv])] <- 0

    pseudospeclist <- list()
    for (i in seq_along(x@pspectra)) {
      # como valor de intensidade de um pico, ele pega o maior valor
      # desse pico entre as amostras
      if (isTRUE(length(x@pspectra[[i]]) > 1) == TRUE) {
        # if database, there is only one value so rowMax doesn't
        # work.
        if (isTRUE(sum(grepl("X", colnames(x@groupInfo)) == TRUE) ==
          0) == TRUE) {
          spc <- cbind(x@groupInfo[x@pspectra[[i]], "mz"], x@groupInfo[
            x@pspectra[[i]],
            sv
          ])
        } else {
          spc <- cbind(x@groupInfo[x@pspectra[[i]], "mz"], rowMax(x@groupInfo[
            x@pspectra[[i]],
            sv
          ]))
        }
        if (isTRUE(sum(is.na(spc)) >= 1) == TRUE) {
          spc <- spc[-(which(is.na(spc), arr.ind = TRUE)[, 1]), ]
        }
      } else {
        spc <- cbind(x@groupInfo[x@pspectra[[i]], "mz"], max(x@groupInfo[
          x@pspectra[[i]],
          sv
        ], na.rm = TRUE))
      }
      pseudospeclist[[i]] <- new("pseudospectrum", id = i, rt = median(x@groupInfo[
        x@pspectra[[i]],
        "rt"
      ]), spectrum = spc)
    }
  }
  ## Possibility to filter pseudospectra by minimum number of peaks
  psvec <- c()
  for (i in seq_along(pseudospeclist)) {
    psvec[i] <- nrow(pseudospeclist[[i]]@spectrum) > min_peaks
  }
  return(pseudospeclist[psvec])
}
