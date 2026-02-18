#' calculateRI
#'
#' This function calculatesRI.
#' Not meant to be called by the user.
#' @param RI Logical (should the software calculate RI), NULL (user will be asked) or path (to the RI file containing 'rt' and 'RI' columns).
#' @param result A list of spectra.
#' @return None
#' @importFrom metaMS addRI write.msp read.msp
#' @importFrom ddpcr quiet
#' @importFrom svDialogs dlg_message
#' @importFrom utils read.csv choose.files
#' @importFrom fritools is_path
#' @export
#' @examples
#' \donttest{
#' result <- metaMS::read.msp((system.file("extdata", "spectra.msp", package = "PipMet")))
#' calculateRI(RI= FALSE, result)
#' }
#' 
calculateRI <- function(RI = NULL, result) {
  if (is.null(RI)) {
    RI <- dlg_message("Add retention index information?", "yesno")$res
    if (RI == "yes") {
      RI <- TRUE
    }
    if (RI == "no") {
      RI <- FALSE
    }
  }
  if (RI == TRUE) {
    again <- "yes"
    while (again == "yes") {
      res <- "yes"
      ri_file <- choose.files()
      if (isEmpty(ri_file) == TRUE & res == "yes") {
        while (res == "yes" & isEmpty(ri_file) == TRUE) {
          RI <- res <- dlg_message(
            "File incorrect or not selected. Add retention index information?",
            "yesno"
          )$res
          if (res == "yes") {
            ri_file <- choose.files()
          } else {
            again <- "no"
          }
        }
      }
      if (isEmpty(ri_file) == FALSE) {
        RI <- ri_file
      }

      if (is_path(RI)) {
        RI <- read.csv(RI)
        if ("rt" %in% colnames(RI) & "RI" %in% colnames(RI)) {
          message("Calculating retention index...")
          quiet(x <- try(result <- addRI(result, RI)))
          if (is(x, "try-error")) {
            again <- dlg_message(
              "File does not meet all requirements for RI calculations. Try again?",
              "yesno"
            )$res
          } else {
            write.msp(result, "spectra_RI.msp", newFile = TRUE)
            dlg_message("The retention index for the spectra was calculated and added to the .msp file.")
            again <- "no"
          }
        } else {
          again <- dlg_message(
            "File must have 'rt' and 'RI' columns. Try again?",
            "yesno"
          )$res
        }
      }
    }
  }
}
