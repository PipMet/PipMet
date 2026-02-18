#' getSpectra
#'
#' This function computes the first subsetp of third step of the workflow. It defines the spectra to proper annotate in NIST MS Search Software.
#' @export
#' @param xdata4 A 'xcmsSet' or 'XCMSnExp' object.
#' @param raw_data A 'XCMSnExp' object.
#' @param min_peaks Numeric. Minimal number of peaks a spectrum must have to be considered a viable spectrum. Default to 5.
#' @param colors A list with colors generated from 'read_data()'.
#' @param column_set Character. Polarity of column used for the chromatography: 'polar', 'non-polar'. If NULL, the user will be asked. Default to NULL.
#' @param prog Character. Configuration of temperature in data acquisition: 'isothermal', 'ramp', 'custom'. If NULL the user will be asked. Default to NULL.
#' @param ion_mode Character. Ion mode acquisition 'positive' or 'negative'. If NULL, the user will be asked. Default to NULL.
#' @param plot_eic Logical. Plot the EIC of each of the 6 most intense m/z in the spectra. Default to FALSE.
#' @return A list with 'xsAnnotate' object with peaks grouped by retention time and correlation peaks information, a 'pseudospectrum' object, a spectra list in .msp format and the ion mode of data acquisition ('negative' or 'positive').
#' @importFrom grDevices dev.off pdf
#' @importFrom stringr str_c
#' @importFrom methods as
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr
#' @importFrom CluMSID writeFeaturelist specplot
#' @importFrom metaMS construct.msp write.msp
#' @importFrom utils memory.limit
#' @importFrom svDialogs dlg_message dlg_input dlg_list


getSpectra <- function(xdata4, raw_data, min_peaks = 5, colors, column_set = NULL,
                       prog = NULL, ion_mode = NULL, plot_eic = FALSE) {
  if (is(xdata4, "XCMSnExp")) {
    xset <- as(xdata4, "xcmsSet")
  } else {
    xset <- xdata4
  }

  if (is.null(ion_mode)) {
    ion_mode <- dlg_list(c("positive", "negative"), multiple = FALSE, title = "Data acquising mode:")$res # ask user for the mode of data acquisition (negative, positive)
  }

  an <- xsAnnotate(xset, polarity = ion_mode)
  anF <- groupFWHM(an, perfwhm = 1)
  rm(an, xset)
  try(memory.limit(1e+05), silent = TRUE)
  anIC <- groupCorr(anF, calcIso = FALSE)
  rm(anF)

  # extract spectra with minimum of 5 peaks
  pslist <- extractSpectra(anIC, min_peaks)

  # create 'pre_anno.csv' where the user annotates the spectra
  # contains id and retention time of each spectra
  writeFeaturelist(pslist)

  # get information about data acquisition
  if (is.null(column_set)) {
    column_set <- dlg_list(c("polar", "non-polar"), multiple = FALSE, title = "Column setup:")$res
  }
  if (is.null(prog)) {
    prog <- dlg_list(c("isothermal", "ramp", "custom"),
      multiple = FALSE,
      title = "Temperatura program:"
    )$res
  }

  # creates a .msp file for spectra
  spectra <- list()
  for (i in seq_len(length(pslist))) {
    x <- data.frame(cbind(pslist[[i]]@spectrum[, 1], (pslist[[i]]@spectrum[
      ,
      2
    ] / max(pslist[[i]]@spectrum[, 2])))) # standardazing intensities by dividing the intensity of each peak from a spectrum by the maximum intensity of that spectra.
    colnames(x) <- c("mz", "into")
    spectra[[i]] <- x
  }
  result <- construct.msp(spectra, extra.info = NULL)
  for (i in seq_len(length(result))) {
    result[[i]]$id <- pslist[[i]]@id
    result[[i]]$rt <- pslist[[i]]@rt
    result[[i]]$Name <- paste0("Unknown ", pslist[[i]]@id)
    result[[i]]$Date <- as.character(Sys.Date())
    result[[i]]$Comments <- paste0("Column class: ", paste0(
      "Standard ",
      column_set
    ), "; ", "ProgramType: ", prog)
    result[[i]]$Ion_mode <- ion_mode
  }
  write.msp(result, "spectra.msp", newFile = TRUE)

  # generate a pdf file with the spectrum, its chromatogram and the XIC
  # of the 6 most intense m/z individually, for every one of the
  # generated spectrum
  if (is.null(plot_eic)) {
    plot_eic <- dlg_message(
      "Plot spectra, with EIC of each of the 6 most intense m/z? (It might take a while)",
      "yesno"
    )$res
  }
  if (plot_eic == "yes" | plot_eic == TRUE) {
    z <- "group"
    rt <- list()
    for (i in seq_len(length(pslist))) {
      rt[[i]] <- as.numeric(pslist[[i]]@rt)
    }
    pdf("EIC_XIC.pdf")
    for (i in seq_len(length(rt))) {
      message(str_c("Printing spectra n.", i, "..."))
      par(mfrow = c(2, 1))
      specplot(pslist[[i]])
      x <- paste0("Unknown ", pslist[[i]]@id)
      crom <- chromatogram(raw_data, rt = c(rt[[i]] - 5, rt[[i]] + 5))
      plot(crom, col = colors[[z]][[2]][colors[[z]][[1]]], main = paste0(
        "Pre-processing - ",
        x
      ))
      legend("right",
        legend = names(colors[[z]][[2]]), col = colors[[z]][[2]],
        fill = colors[[z]][[2]], box.lty = 0, cex = 0.8, bg = "transparent"
      )
      # plot EIC from the 6 most intense ion-fragm
      par(mfrow = c(2, 3))
      f <- pslist[[i]]@spectrum[order(pslist[[i]]@spectrum[, 2], decreasing = TRUE)]
      # sort(pslist[[i]]@spectrum[, 1], decreasing = TRUE)
      for (ii in seq_len(6)) {
        crom <- chromatogram(raw_data, rt = c(rt[[i]] - 5, rt[[i]] +
          5), mz = c(as.numeric(f[[ii]]) - 0.6, as.numeric(f[[ii]]) +
          0.6))
        plot(crom, col = colors[[z]][[2]][colors[[z]][[1]]], main = paste0(
          "EIC - mz",
          f[[ii]]
        ))
        legend("right",
          legend = names(colors[[z]][[2]]), col = colors[[z]][[2]],
          fill = colors[[z]][[2]], box.lty = 0, cex = 0.8, bg = "transparent"
        )
      }
    }
    dev.off()
  }


  return(list(
    anIC = anIC, pslist = pslist, result = result, ion_mode = ion_mode,
    column_set = column_set, prog = prog
  ))
}
