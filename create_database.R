#' In-house database from processed files
#'
#' Function for creating a in-house database from processed GC-MS experiment files. The user may choose to add more information, such as Formula, Monoisotopic Mass, CAS, ChemSpider, InChIKey, PubChem ID, Class and retention index. The function is capable of search and retrieving most of the information from the online databases such as PubChem and NIST uatomatically through 'webchem' package or add manually the information.
#' @keywords database spectra
#' @export
#' @param apslist List of spectra with annotation. Only annotated spectra will be use for the database construction.
#' @param ion_mode Character. Ion mode of data acquisition ('negative' or 'positive').
#' @param column_set Character. Polarity of column used for the chromatography: 'polar', 'non-polar'. If NULL, the user will be asked. Default to NULL.
#' @param prog Character. Configuration of temperature in data acquisition: 'isothermal', 'ramp', 'custom'. If NULL the user will be asked. Default to NULL.
#' @param ion_mode Character. Ion mode acquisition 'positive' or 'negative'. If NULL, the user will be asked. Default to NULL.
#' @param info Logical. If TRUE, the user cann add more information to the creatin of the database, such as CAS number, PubMed and ChemSpider ID and InChiKey, among others. If NULL, the user will be asked. Default to NULL.
#' @param Ri Character. Retention index calculation method. Supported are 'lee', 'linear', 'kovats' and 'alcane'. If NULL, the user will be asked. Default ot NULL.
#' @return None.
#' @importFrom svDialogs dlg_message dlg_list
#' @importFrom utils read.csv write.csv
#' @importFrom metaMS write.msp
#' @importFrom webchem get_cid pc_prop cts_convert nist_ri
#' @importFrom methods .hasSlot
#' @importFrom pracma isempty
#' @importFrom ddpcr quiet
#' @importFrom stringr str_to_lower
#' @examples
#' \donttest{
#' \dontrun{
#' load(system.file('extdata', 'apslist.RData', package = 'PipMet'))
#' create_database(apslist, 
#'                  column_set = 'non-polar', 
#'                  prog = 'ramp', 
#'                  ion_mode = 'positive', 
#'                  info = FALSE,
#'                  Ri = 'kovats')
#' }
#' }
create_database <- function(apslist, column_set = NULL, prog = NULL, ion_mode = NULL,
    info = NULL, Ri = NULL) {

    # get only annotated spectra
    pslist <- list()
    count <- 1
    names <- vector()
    for (i in seq_along(apslist)) {
        if (!isempty(apslist[[i]]@annotation)) {
            quiet(pslist[count] <- apslist[[i]])
            count <- count + 1
            names <- c(names, apslist[[i]]@annotation)
        }
    }

    # add information from internet?
    if (is.null(info)) {
        info <- dlg_message("Would you like to add more informations about the compounds?",
            type = "yesno")$res
    }
    if (info == "yes" | info == TRUE) {
        dataInfo <- matrix(nrow = length(pslist), ncol = 10)
        colnames(dataInfo) <- c("Name", "formula", "exact.mass", "rt", "CAS",
            "ChemSpider", "InChIKey", "PubChem ID", "Class", "RI")
        dataInfo <- as.data.frame(dataInfo)
        for (i in seq_len(length(pslist))) {
            dataInfo[i, 1] <- pslist[[i]]@annotation
            dataInfo[i, "rt"] <- pslist[[i]]@rt
        }
        write.csv(dataInfo, "dataInfo.csv", row.names = FALSE, na = "")
        dlg_message("A \"dataInfo.csv\" file was written in the project folder. Please fill out the fields and press \"ok\".",
            "ok")$res
        dataInfo <- read.csv("dataInfo.csv", na = "", check.names = FALSE)
        for (i in seq_len(length(pslist))) {
            if (.hasSlot(pslist[[i]], "RI") == TRUE) {
                dataInfo[i, "RI"] <- pslist[[i]]@RI
            }
            dataInfo[i, "PubChem ID"] <- get_cid(dataInfo[i, 1])$cid
            f <- pc_prop(dataInfo[i, "PubChem ID"])
            if (!is.na(f["ChemSpiderID"])) {
                dataInfo[i, "formula"] <- f[1, "MolecularFormula"]
                dataInfo[i, "InChIKey"] <- f[1, "InChIKey"]
                dataInfo[i, "exact.mass"] <- f[1, "MonoisotopicMass"]
            }
        }
        dataInfo[, "CAS"] <- unlist(cts_convert(dataInfo[, "PubChem ID"], "PubChem CID",
            "cas", match = "first"))
        # ask information for retention index search
        if (dlg_message("Look for retention index in NIST?", type = "yesno")$res ==
            "yes") {
            if (is.null(Ri)) {
                Ri <- dlg_list(c("kovats", "linear", "alkane", "lee"), multiple = FALSE,
                  title = "Retention time index:")$res
            }
            if (is.null(column_set)) {
                column_set <- dlg_list(c("polar", "non-polar"), multiple = FALSE,
                  title = "Column setup:")$res
            }
            if (is.null(prog)) {
                prog <- dlg_list(c("isothermal", "ramp", "custom"), multiple = FALSE,
                  title = "Temperature program:")$res
            }
            x <- dlg_list(c("CAS", "Name", "InChIKey"), multiple = FALSE, title = "Search for retention index based on:")$res
            for (i in seq_len(length(pslist))) {
                if (!is.na(dataInfo[i, x])) {
                  dataInfo[i, "RI"] <- mean(nist_ri(dataInfo[i, x], from = str_to_lower(x),
                    type = Ri, polarity = column_set, temp_prog = prog)$RI)
                } else {
                  dataInfo[i, "RI"] <- mean(nist_ri(dataInfo[i, "Name"], from = "name",
                    type = Ri, polarity = column_set, temp_prog = prog)$RI)
                }
            }
        }
        # write informations for checking
        write.csv(dataInfo, file = "DatabaseInfo.csv", row.names = FALSE)
        dlg_message("Check and fill the 'DatabaseInfo.csv' file in myDir and press 'OK'.")$res
        dataInfo <- read.csv("DatabaseInfo.csv", na.string = c("NA", ""), colClasses = "character",
            sep = ",")
        if (!sum(grep("NA", dataInfo$Name)) == 0) {
            dataInfo <- dataInfo[-which(is.na(dataInfo$Name)), ]
        }
    }

    # create spectra in .msp file format
    spectra <- list()
    for (i in seq_len(length(pslist))) {
        x <- cbind(pslist[[i]]@spectrum[, 1], (pslist[[i]]@spectrum[, 2]/max(pslist[[i]]@spectrum[,
            2])))  # pslistronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
        x <- data.frame(x)
        colnames(x) <- c("mz", "into")
        spectra[[i]] <- x
        rm(x)
    }
    result <- metaMS::construct.msp(spectra, extra.info = NULL)
    for (i in seq_along(result)) {
        if (info == "yes") {
            result[[i]]$rt <- dataInfo[i, "rt"]
            result[[i]]$Name <- dataInfo[i, "Name"]
            result[[i]]$Formula <- dataInfo[i, "formula"]
            result[[i]]$MW <- dataInfo[i, "exact.mass"]
            result[[i]]$CAS <- dataInfo[i, "CAS"]
            result[[i]]$ChemSpiderID <- dataInfo[i, "ChemSpiderID"]
            result[[i]]$InChIKey <- dataInfo[i, "InChIKey"]
            result[[i]]$PubChemID <- dataInfo[i, "PubChem ID"]
            result[[i]]$Class <- dataInfo[i, "Class"]
            result[[i]]$RI <- dataInfo[i, "RI"]
            result[[i]]$Ion_mode <- ion_mode
            result[[i]]$Comments <- paste0("Column class: ", paste0("Standard ",
                column_set), "; ", "ProgramType: ", prog, "; RI ", Ri, " :",
                result[[i]]$RI <- dataInfo[i, "RI"])
        }
        if (.hasSlot(pslist[[1]], "RI") == TRUE) {
            result[[i]]$RI <- pslist[[i]]@RI
        }
        result[[i]]$Date <- as.character(Sys.Date())
        result[[i]]$Comments <- paste0("Column class: ", paste0("Standard ",
            column_set), "; ", "ProgramType: ", prog)
    }
    names(result) <- names
    metaMS::write.msp(result, file = paste0("Database_", Sys.Date(), ".msp"),
        newFile = TRUE)
    dlg_message("Database creation done!")$res
}
