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
#'    pictures = FALSE,
#'    example = TRUE
#' )
#' }}
read_data <- function(peakMonitor = NULL, ions = NULL, myDir = NULL, sample_dir = NULL,
                      metadata = NULL, extension = NULL, pictures = TRUE, example = FALSE, 
                      pic_extension = c(".tiff", ".png"), sample_names = NULL) {

  if (example == FALSE) {
    # --- MODO NORMAL ---
    if (is.null(sample_dir) | missing(sample_dir)) {
      sample_dir <- choose.dir(default = getwd(), caption = "Please, select the Samples directory")
    }
    sample_dir <- normalizePath(sample_dir, winslash = "/", mustWork = FALSE)

    if (is.null(myDir) | missing(myDir)) {
      myDir <- dlgInput("Name your project", Sys.info()["user"])$res
    }
    # Cria a pasta no diretório atual
    if (!dir.exists(myDir)) dir.create(myDir, showWarnings = FALSE)
    myDir <- normalizePath(myDir, winslash = "/", mustWork = FALSE)

    if (is.null(metadata)) {
      if (is.null(extension)) {
        extension <- dlg_list(c(".mzML", ".mzXML"), multiple = FALSE, title = "Files extension:")$res
      }
      
      metd <- dlg_message("Metadata table already exists?", "yesno")$res
      if (metd == "yes") {
        meta_file <- choose.files(caption = "Select Metadata File", multi = FALSE)
        metadata <- read.csv(meta_file, na.string = c("NA", ""), colClasses = "character")
        metadata$file <- sapply(metadata$file, function(x) file.path(sample_dir, basename(x)))
      } else {
        files <- list.files(sample_dir, full.names = FALSE, pattern = extension, recursive = TRUE)
        metadata <- data.frame(sample = sub(basename(files), pattern = extension, replacement = "", fixed = TRUE),
                               group = "A", class = "Sample", tec_rep = "1", bio_rep = "1", 
                               file = file.path(sample_dir, files), stringsAsFactors = FALSE)
        write.csv(metadata, file.path(myDir, "metadata.csv"), row.names = FALSE)
        dlg_message("Fill 'metadata.csv' in your project folder and press OK.", type = "ok")
      }
    }
  } else {
    # --- MODO EXEMPLO ---
    message("Configuring example mode...")
    meta_path <- system.file("extdata", "metadata.csv", package = "PipMet")
    metadata <- read.csv(meta_path, na.string = c("NA", ""), colClasses = "character")
    
    # Ajusta extensões e caminhos internos do pacote
    metadata$file <- sapply(metadata$file, function(x) {
      fname <- basename(x)
      if (!grepl("\\.mzXML$|\\.mzML$", fname, ignore.case = TRUE)) fname <- paste0(fname, ".mzXML")
      system.file("extdata", fname, package = "PipMet")
    })

    # Pergunta o nome da pasta de exemplo para o usuário
    proj_name <- dlgInput("Name your example project folder", "PipMet_Example_Test")$res
    myDir <- file.path(getwd(), proj_name)
    if (!dir.exists(myDir)) dir.create(myDir, showWarnings = FALSE)
    myDir <- normalizePath(myDir, winslash = "/", mustWork = FALSE)
    
    message("Project folder created at: ", myDir)
  }

  # --- PROCESSAMENTO COMUM ---
  
  # Cores
  if (nrow(metadata) > 1) {
    x_cols <- colnames(metadata)[!colnames(metadata) %in% c("file", "all", "all2")]
    colors <- vector(mode = "list", length = length(x_cols))
    names(colors) <- x_cols
    for (i in seq_along(x_cols)) {
      levs <- unique(metadata[[x_cols[i]]])
      pal <- if(length(levs) <= 9) RColorBrewer::brewer.pal(max(3, length(levs)), "Set1")[1:length(levs)] else rainbow(length(levs))
      names(pal) <- levs
      colors[[i]] <- list(metadata[[x_cols[i]]], pal)
      names(colors[[i]]) <- c(x_cols[i], paste0(x_cols[i], "_colors"))
    }
  } else { colors <- list() }

  # Leitura XCMS
  message("Reading MS files... this may take a minute.")
  raw_data <- readMSData(metadata$file, pdata = new("AnnotatedDataFrame", metadata), mode = "onDisk")

  # Visualização
  if (pictures) {
    message("Generating visualization plots...")
    vis_dir <- file.path(myDir, "Visualization_results")
    if (!dir.exists(vis_dir)) dir.create(vis_dir)
    
    # Salva diretório atual e garante retorno
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(vis_dir)

    bpc <- chromatogram(raw_data, aggregationFun = "max")
    tic <- chromatogram(raw_data, aggregationFun = "sum")

    for (i in seq_along(colors)) {
      group_name <- names(colors)[i]
      sample_cols <- colors[[i]][[2]][metadata[[group_name]]]
      
      if (".png" %in% pic_extension) {
        png(paste0(group_name, "_chromatograms.png"), width = 1000, height = 1000, res = 120)
        par(mfrow = c(2, 1))
        plot(bpc, col = sample_cols, main = paste("BPC by", group_name))
        plot(tic, col = sample_cols, main = paste("TIC by", group_name))
        dev.off()
      }
    }
    message("Plots saved in: ", vis_dir)
  }

  setwd(myDir)
  return(list(colors = colors, metadata = metadata, raw_data = raw_data, myDir = myDir))
}