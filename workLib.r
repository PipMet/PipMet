#' workLib
#'
#' Main function for GC-MS standards processing for internal library development.
#' @keywords standards, library, internal
#' @export
#' @return None.
#' @param myDir Path to working directory. Default to none.
#' @param libname Character. Name of the library. Default to NULL. If NULL the user will be asked.
#' @param sample_dir Path to sample directory. Default to none.
#' @param lib_metadata Path to .csv file containing metadata. 'compound', 'formula', 'exact.mass', 'rt', 'file', 'CAS', 'ChemSpider', 'class', 'RI', 'InChIKey' columns must be included. Default to none.
#' @param extension Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'. Default to none.
#' @param example Logical. If is example, pop-ups won't appear.
#' @param parallel Character. Sort of parallelization for code to perfom. Supported are 'Serial Param', 'Snow Param', 'MultiCore Param'. For more information, check the BiocParallel R package. Default to NULL. If parallel = NULL, the user will be asked.
#' @param column_set Character. Polarity of column used for the chromatography: 'polar', 'non-polar'. If NULL, the user will be asked. Default to NULL.
#' @param prog Character. Configuration of temperature in data acquisition: 'isothermal', 'ramp', 'custom'. If NULL the user will be asked. Default to NULL.
#' @param ion_mode Character. Ion mode acquisition 'positive' or 'negative'. If NULL, the user will be asked. Default to NULL.
#' @param instrument_type Character. Instrument of acquisition. Default to NULL. If NULL, the user will be asked.
#' @param cores Numeric. Number of cores to be used in Snow Param. Default to NULL. If NULL, the user will be asked. Set cores = 0 to Serial Param.
#' @param Ri_info Character. If the user wants to add retention information or not. Supported are 'From file', 'Retrieve from NIST', 'No RI information'. Default to NULL
#' @param Ri Character. Retention index calculation method. Supported are 'lee', 'linear', 'kovats' and 'alcane'. If NULL, the user will be asked. Default ot NULL.
#' @param RI Path to retention index information .csv file. Only if Ri_info == 'From file'. Default to NULL. If NULL the user will be asked.
#' @importFrom methods as
#' @importFrom ddpcr quiet
#' @importFrom stringr str_c
#' @importFrom utils choose.dir menu read.csv write.csv
#' @importFrom metaMS addRI write.msp construct.msp
#' @importFrom MSnbase readMSData
#' @importFrom ProtGenerics filterRt
#' @importFrom tools file_path_as_absolute
#' @importFrom xcms MatchedFilterParam findChromPeaks
#' @importFrom BiocParallel register SerialParam
#' @importFrom webchem nist_ri as.cas
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr
#' @importFrom tcltk tkmessageBox
#' @importFrom ddpcr quiet
#' @importFrom fritools is_path
#' @importFrom svDialogs dlgInput dlg_message
#' @examples
#' \donttest{
#' \dontrun{
#' workLib(
#'   sample_dir = system.file('extdata', package = 'PipMet'),
#'   lib_metadata = system.file('extdata', 'lib_metadata.csv', package = 'PipMet'),
#'   extension = '.mzXML',
#'   myDir = '~/',
#'   example = TRUE
#' )
#' }
#' }
workLib <- function(myDir = NULL, libname = NULL, sample_dir = NULL, lib_metadata = NULL,
    extension = NULL, example = FALSE, parallel = NULL, ion_mode = NULL, prog = NULL,
    column_set = NULL, instrument_type = NULL, cores = 1, Ri = NULL, RI = NULL,
    Ri_info = NULL) {

  # parallel setup

  setupParallel <- function(parallel, cores) {
    if (is.null(parallel)) {
      parallel <- svDialogs::dlg_list(c("Serial Param", "Snow Param", "MultiCore Param"),
                                     multiple = FALSE, title = "Parallelization mode:")$res
    }
    if (parallel == "Serial Param") {BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)}
    if (parallel == "Snow Param") {
      if (is.null(cores)) {
        cores <- as.numeric(svDialogs::dlgInput(paste0("Number of cores used (you have ", 
                                                      parallel::detectCores(), " cores available):"), 
                                               parallel::detectCores() - 2)$res)
      }
      BiocParallel::register(BiocParallel::SnowParam(workers = cores), default = TRUE)
    }
    if (parallel == "MultiCore Param") BiocParallel::register(BiocParallel::MulticoreParam(), default = TRUE)
    return(parallel)
  }
  
  parallel <- setupParallel(parallel, cores)


    # get directory for samples and project and set library metadata
    if (!example == TRUE) {

        # lib name
        if (is.null(libname)) {
            libname <- dlgInput("Name your library", "Library X")$res
        }

        # samples directory
        if (is.null(sample_dir)) {
            sample_dir <- choose.dir(default = getwd(), caption = "Please, select the Samples directory, should be C:/Users/_/Samples")
        }

        # project directory
        if (is.null(myDir)) {
            myDir <- dlgInput("Name your project", Sys.info()["user"])$res
            if (!dir.exists(myDir) == TRUE) {
                dir.create(myDir)
            }
        }
        setwd(myDir)
        myDir <- getwd()

        # files extension
        if (is.null(extension)) {
            extension <- dlg_list(c(".mzML", ".mzXML"), multiple = TRUE, title = "Files extension:")$res
        }

        # library metadata set up
        if (is.null(lib_metadata)) {
            metd <- dlg_message("Metadata table already exists?", "yesno")$res
            if (metd == "yes") {
                # metadata already exists
                lib_metadata <- read.csv(choose.files(), na.string = c("NA",
                  ""), colClasses = "character", sep = ",")
                if (!"file" %in% colnames(lib_metadata)) {
                  lib_metadata$file <- paste0(sample_dir, "/", lib_metadata$compound,
                    extension)
                } else {
                  for (i in seq_len(nrow(lib_metadata))) {
                    # check if they are just filenames or filepath (they
                    # need to be path)
                    if (!is_path(lib_metadata$file[i])) {
                      lib_metadata$file[i] <- paste0(sample_dir, "/", basename(lib_metadata$file[i]))
                    }
                  }
                }
            } else {
                # metadata doesn't exist
                files <- list.files(sample_dir, full.names = TRUE, pattern = extension,
                  recursive = TRUE)
                lib_metadata <- matrix(nrow = length(files), ncol = 10)
                colnames(lib_metadata) <- c("compound", "formula", "exact.mass",
                  "rt", "file", "CAS", "ChemSpider", "class", "RI", "InChIKey")
                lib_metadata[, "file"] <- files
                lib_metadata[, "compound"] <- sub(basename(files), pattern = extension,
                  replacement = "", fixed = TRUE)
                write.csv(lib_metadata, "lib_metadata.csv", row.names = FALSE)
                dlg_message("A file 'lib_metadata.csv' was created in your directory. Fill the sheet before continuing. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.",
                  "yesno")
                while (file.exists("lib_metadata.csv") == FALSE) {
                  write.csv(lib_metadata, "lib_metadata.csv", row.names = FALSE)
                  dlg_message("A file 'lib_metadata.csv' was created in you directory. Fill the sheet before continuing. You can create new columns to describe samples, such as 'strain'. After filling the sheet, press 'ok'.",
                    type = "ok")
                  lib_metadata <- read.csv("lib_metadata.csv", na.string = c("NA",
                    ""), colClasses = "character", sep = ",", dec = ".")
                }
                lib_metadata <- read.csv("lib_metadata.csv", na.string = c("NA",
                  ""), colClasses = "character", sep = ",", dec = ".")
            }
        }
        if (!is.null(lib_metadata)) {
            # if lib_metadata is provided if metadata if a path to
            # lib_metadata file
            if (is_path(lib_metadata)) {
                lib_metadata <- read.csv(lib_metadata, na.string = c("NA", ""),
                  colClasses = "character", sep = ",")
                if (!"file" %in% colnames(lib_metadata)) {
                  lib_metadata$file <- paste0(sample_dir, "/", lib_metadata$compound,
                    extension)
                } else {
                  for (i in seq_len(nrow(lib_metadata))) {
                    # check if they are just filenames or filepath (they
                    # need to be path)
                    if (!is_path(lib_metadata$file[i])) {
                      lib_metadata$file[i] <- paste0(sample_dir, "/", basename(lib_metadata$file[i]))
                    }
                  }
                }
            }
            while (!is(lib_metadata, "data.frame") | !is(lib_metadata, "matrix")) {
                # lib_metadata may also be a data.frame or matrix. if none,
                # ask for lib_metadata again
                dlg_message("Library metadata format not supported. Upload a .csv file with, at least, 'compound', 'rt' and 'file' columns!")$res
                lib_metadata <- read.csv(choose.files(), na.string = c("NA",
                  ""), colClasses = "character", sep = ",")
            }
        }
        # if example = TRUE
    } else {
        sample_dir <- system.file("extdata", package = "PipMet")
        lib_metadata <- read.csv(system.file("extdata", "lib_metadata.csv", package = "PipMet"),
            na.string = c("NA", ""), colClasses = "character", sep = ",", dec = ".")
        extension <- ".mzML"
        for (i in seq_len(nrow(lib_metadata))) {
            lib_metadata$file[i] <- system.file("extdata", lib_metadata$file[i],
                package = "PipMet")
        }
        myDir <- dlgInput("Name your project", "PipMet_example")$res
        lib_name <- "Metabolite_mix"
        setwd("~/")
        dir.create(myDir)
        setwd(myDir)
    }

    # remove empty lines/columns from lib_metadata
    lib_metadata <- Filter(function(x) !all(is.na(x)), lib_metadata)  # remove empty columns
    c <- vector()
    for (i in seq_len(nrow(lib_metadata))) {
        if (is.na(lib_metadata$compound[i])) {
            c <- c(c, i)
        }
    }
    if (!is(c, "logical")) {
        lib_metadata <- lib_metadata[-c, ]
    }


    # begin processing

    # get mode of data acquisiton
    if (is.null(ion_mode)) {
        ion_mode <- dlg_list(c("positive", "negative"), multiple = FALSE, title = "Ion mode of data acquisition:")$res
    }

    # read and filter files to a 5s window
    message(str_c("Reading ", length(unique(lib_metadata$file)), " files..."))
    raw_data <- list()
    for (i in seq_len(nrow(lib_metadata))) {
        quiet(raw_data[[i]] <- readMSData(lib_metadata$file[i], mode = "onDisk"))
        quiet(raw_data[[i]] <- filterRt(raw_data[[i]], c(as.numeric(lib_metadata[i,
            "rt"]) - 5, as.numeric(lib_metadata[i, "rt"]) + 5)))
    }

    # peak picking
    message(str_c("Preprocessing ", length(unique(lib_metadata$file)), " files..."))
    quiet(xdata <- lapply(raw_data, FUN = findChromPeaks, param = MatchedFilterParam(fwhm = 5,
        binSize = 0.5, steps = 2, mzdiff = 0.5, snthresh = 2, max = 500)))

    message("Grouping peaks into spectra...")
    # spectra defition
    quiet(xset <- lapply(xdata, FUN = as, Class = "xcmsSet"))

    # create a xsAnnotate object
    quiet(an <- lapply(xset, FUN = xsAnnotate, polarity = ion_mode))

    # group by rt
    quiet(anF <- lapply(an, FUN = groupFWHM, perfwhm = 1))

    # group by correlation information
    quiet(anIC <- lapply(anF, FUN = groupCorr, calcIso = FALSE))

    # extract to spectra list
    quiet(pslist <- lapply(anIC, FUN = extractSpectra, min_peaks = 5))

    names(pslist) <- names(anIC) <- names(an) <- names(anF) <- names(xset) <- names(xdata) <- names(raw_data) <- lib_metadata[,
        "compound"]

    # get information on retention index, column and temperature program
    if (is.null(Ri_info)) {
        Ri_info <- dlg_list(c("From file", "Retrieve from NIST", "No RI information"),
            multiple = FALSE, title = "Retention time index: ")$res
    }
    if (!Ri_info == "No Ri information") {
        if (is.null(Ri)) {
            Ri <- dlg_list(c("kovats", "linear", "alkane", "lee"), multiple = FALSE,
                title = "Retention time index")$res
        }
    }
    if (Ri_info == "From file") {
        if (is.null(RI)) {
            RI <- read.csv(choose.files())
        }
        if (is_path(RI)) {
            RI <- read.csv(RI)
        }
    }

    if (is.null(column_set)) {
        column_set <- dlg_list(c("polar", "non-polar"), multiple = FALSE)$res
    }
    if (is.null(prog)) {
        prog <- dlg_list(c("isothermal", "ramp", "custom"), multiple = FALSE)$res
    }
    if (is.null(instrument_type)) {
        instrument_type <- dlg_input("Type of instrument of acquisition", "GC-EI-Q")$res
    }

    if (Ri_info == "Retrieve from NIST") {
        message("Retrieving RI from CAS numbers, if available...")
        # retrieve RI from CAS numbers (from NIST database), if RI is not
        # present
        lib_metadata[, "RI_"] <- NA
        for (i in seq_len(nrow(lib_metadata))) {
            lib_metadata[i, "CAS"] <- gsub(pattern = "-", replacement = "", lib_metadata[i,
                "CAS"])  # remove hifens from CAS
            if (is.na(lib_metadata[i, "RI_"])) {
                quiet(lib_metadata[i, "RI_"] <- mean(nist_ri(as.cas(lib_metadata[i,
                  "CAS"]), from = "cas", type = Ri, polarity = column_set, temp_prog = prog)$RI))
            }
        }
    }

    # create .msp structure and files
    message(str_c("Writing ", length(pslist), " spectra for manual validation..."))

    if (!dir.exists("spectra")) {
        dir.create("spectra")
    }
    setwd("spectra")

    for (i in seq_len(length(pslist))) {
        spectra <- list()
        for (ii in seq_len(length(pslist))) {
            x <- cbind(pslist[[i]][[ii]]@spectrum[, 1], (pslist[[i]][[ii]]@spectrum[,
                2]/max(pslist[[i]][[ii]]@spectrum[, 2])))  # padronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
            x <- data.frame(x)
            colnames(x) <- c("mz", "into")
            spectra[[ii]] <- x
        }
        result <- construct.msp(spectra, extra.info = NULL)
        for (ii in seq_len(length(result))) {
            result[[ii]]$Name <- paste0(pslist[[i]][[ii]]@id, " - Candidate to ",
                names(pslist)[i])
            result[[ii]]$id <- pslist[[i]][[ii]]@id
            result[[ii]]$rt <- pslist[[i]][[ii]]@rt
            result[[ii]]$RI <- lib_metadata$RI[[i]]
        }
        metaMS::write.msp(result, paste0(names(pslist)[i], "_", Sys.Date(), ".msp"),
            newFile = TRUE)
    }

    # alerts the user that a .msp file for each component of the
    # lib_metadata list was created, so the user uploads once at a time to
    # NIST MSSearch, validate and answer the folowwing pop-up asking which
    # index of mass spectra is the component tcltk::tkmessageBox(title =
    # 'Library development', message = 'For each file processed, a '.msp'
    # file was created in the files directory. Upload to NIST MS Search and
    # select the mass spectrum correspondent to the standard.', icon =
    # 'info', type = 'ok')
    dlg_message("For each file processed, a '.msp' file was created in the files directory. Upload to NIST MS Search and select the mass spectrum correspondent to the standard.")$res

    # ppslist is the pslist corrected
    ppslist <- pslist

    # validating spectra

    done <- 2
    while (done == 2) {
        for (i in seq_len(length(pslist))) {
            id <- dlg_list(c(seq(length(pslist[[i]])), "None"), multiple = FALSE,
                title = names(pslist)[[i]])$res
            if (id == "None") {
                ppslist[[i]] <- "NA"
            } else {
                ppslist[[i]] <- pslist[[i]][[as.numeric(id)]]
            }
        }
        done <- menu(c("Continue", "Repeat"), graphics = TRUE, title = "Validation done!")
    }

    message(str_c("Writing ", length(ppslist), " validated spectra..."))
    # generate final .msp file
    spectra <- list()
    for (i in seq_len(length(ppslist))) {
        x <- cbind(ppslist[[i]]@spectrum[, 1], (ppslist[[i]]@spectrum[, 2]/max(ppslist[[i]]@spectrum[,
            2])))  # padronizo dividindo todas as intensidades de um mesmo espectro pela maior intensidade no mesmo (fica tipo 1 e 0,X ou seja, porcentagens)
        x <- data.frame(x)
        colnames(x) <- c("mz", "into")
        spectra[[i]] <- x
        rm(x)
    }
    result <- metaMS::construct.msp(spectra, extra.info = NULL)
    lib_metadata[is.na(lib_metadata)] <- ""
    for (i in seq_len(length(result))) {
        result[[i]]$id <- ppslist[[i]]@id
        result[[i]]$rt <- ppslist[[i]]@rt
        result[[i]]$Name <- names(ppslist)[i]
        result[[i]]$Formula <- lib_metadata[grep(names(ppslist)[i], lib_metadata[,
            1], value = FALSE), "formula"]
        result[[i]]$MW <- lib_metadata[grep(names(ppslist)[i], lib_metadata[,
            1], value = FALSE), "exact.mass"]
        result[[i]]$CAS <- lib_metadata[grep(names(ppslist)[i], lib_metadata[,
            1], value = FALSE), "CAS"]
        result[[i]]$Class <- lib_metadata[grep(names(ppslist)[i], lib_metadata[,
            1], value = FALSE), "class"]
        result[[i]]$Date <- as.character(Sys.Date())
        result[[i]]$Instrument_type <- instrument_type
        result[[i]]$Comments <- paste0("Column class: ", paste0("Standard ",
            column_set), "; ", "ProgramType: ", prog, "; ", "RI: ", Ri)
        result[[i]]$Ion_mode <- ion_mode
        if ("RI" %in% colnames(lib_metadata)) {
            result[[i]]$RI <- lib_metadata[grep(names(ppslist)[i], lib_metadata[,
                1], value = FALSE), "RI"]
        }
        if ("InChIKey" %in% colnames(lib_metadata)) {
            result[[i]]$InChIKey <- lib_metadata[grep(names(ppslist)[i], lib_metadata[,
                1], value = FALSE), "InChIKey"]
        }
    }
    names(result) <- names(ppslist)
    setwd(myDir)

    if (Ri_info == "From file") {
        message("Calculating retention index...")
        result <- addRI(result, RI)
    }

    metaMS::write.msp(result, paste0("Library_", libname, ".msp"), newFile = TRUE)

    dlg_message("Done! Internal library development finalized. A '.msp' file was created in your directory for conversion to NIST MS Search Library.")
}


