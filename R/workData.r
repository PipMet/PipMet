#' workData
#'
#' Main function for GC data processing. It read files in '.mzXML' and '.mzML' formats, apply algorithms for for peacking, grouping, retention time correction and filling missing peaks. The next step proposes spectra based on retention time window and correlation information. After annotation, the function performs quantification by normalizing data and apply t test and principal components analysis, plotting pictures for the evaluation of the steps performed. At the end, the user can choose to generate a internal library with the identified compounds to be uploaded into NIST MS Search software.
#' @keywords metadata in-house database preprocessing spectra spectrum peakpicking grouping RI rt
#' @export
#' @param myDir Path to working directory. Default to none.
#' @param sample_dir Path to sample directory. Default to none.
#' @param metadata Path to .csv file or an R data.frame object containing metadata. At least 'sample' and 'file' columns must be included. Default to none.
#' @param extension Extension of mass spectrometry files to read. Only accepted '.mzML' and '.mzXML'. Default to none.
#' @param pictures Logical. If pictures should be plotted or not.
#' @param example Logical. If example = TRUE, the metadata and other needed files will be loaded from package files.
#' @param filter Numeric. Intensity threshold for the peak detection. Default to NULL. When NULL, the user will be asked for a number. Set filter = 0 for no intensity filtering.
#' @param peakMonitor Logical. Are there peak to monitor throuhout the workflow? If NULL, the user will be asked. Default to NULL.
#' @param Ri Character. Retention index calculation method. Supported are 'lee', 'linear', 'kovats' and 'alcane'. If NULL, the user will be asked. Default ot NULL.
#' @param parallel Character. Sort of parallelization for code to perfom. Supported are 'Serial Param', 'Snow Param', 'MultiCore Param'. For more information, check the BiocParallel R package. Default to NULL. If parallel = NULL, the user will be asked.
#' @param pic_extension Character. Pictures format to generate. Supported = '.tiff', '.png'. Default to c('.tiff', '.png').
#' @param group Character. Name from 'metadata' column names to group the samples. Default to 'group'.
#' @param derivatization Character. Kind of derivatization the samples were prepared with. Supported are 'Trimethylsilyl' and 'None'. If NULL, the user will be asked. Default to 'NULL'.
#' @param cores Numeric. Number of cores to be used in Snow Param. Default to NULL. If NULL, the user will be asked. Set cores = 0 to Serial Param.
#' @param column_set Character. Polarity of column used for the chromatography: 'polar', 'non-polar'. If NULL, the user will be asked. Default to NULL.
#' @param prog Character. Configuration of temperature in data acquisition: 'isothermal', 'ramp', 'custom'. If NULL the user will be asked. Default to NULL.
#' @param RI Logical or path to retention index .csv file. Addition of retention index information to the spectra. If RI = TRUE, the user will be asked to provide a .csv file with 'rt' and 'RI' columns. If RI = path to the .csv files, the retention index will be calculated. Default to FALSE.
#' @param ion_mode Character. Ion mode acquisition 'positive' or 'negative'. If NULL, the user will be asked. Default to NULL.
#' @param plot_eic Logical. Plot the EIC of each of the 6 most intense m/z in the spectra. Default to NULL.
#' @param replicate Character or Logical. If FALSE, there is no replicate informations in the metadata table. Otherwise, inform the name in the metadata column containing replicate information. If NULL, the user will be asked. Default to NULL.
#' @param removeCompounds Logical. If TRUE, the user may choose from a pop-up identified compounds to remove from the quantification table. If NULL, the user will be asked. Default to 'NULL'.
#' @param mergeCompounds Logical. If TRUE, the user may choose from a pop-up what compounds identified should be representated as one with intensities summed. Exemple: derivatizations derivatives. If NULL, the user will be asked. Default to 'NULL'.
#' @param pre_anno Path to pre_anno.csv annotation file or pre_anno table. If NULL, the pre_anno.csv in the folder will be read. Default to NULL.
#' @param min_peaks Numeric. Minimal number of peaks a spectrum must have to be considered a viable spectrum. Default to 5.
#' @param sample_names Character. Name of metadata column to name the samples. Default to NULL. If NULL, the user will be asked.
#' @return A list containing (1) the path of working folder, (2) the metadata table, (3) the annotated pseudospectra list, (4) a OnDiskMSnExp object, (5) a XCMSnExp or xcmsSet object, (6) a xsAnnotate object, (7) a list of colors used and (8) the normalized instensities quantification table
#' @importFrom methods as new
#' @importFrom svDialogs dlgInput dlg_message dlg_list dlg_confirm
#' @importFrom ddpcr quiet
#' @importFrom utils choose.dir menu read.csv write.csv write.table choose.files
#' @importFrom metaMS addRI write.msp
#' @importFrom parallel detectCores
#' @importFrom BiocParallel register SerialParam SnowParam MulticoreParam
#' @importFrom tcltk tkmessageBox
#' @importFrom fritools is_path
#' @examples
#' \donttest{
#' \dontrun{
#' result <- workData(
#'    extension = '.mzXML',
#'    myDir = '~/',
#'    example = TRUE,
#'    pictures = TRUE
#' )
#' }
#' }
workData <- function(myDir = NULL, sample_dir = NULL, metadata = NULL, extension = NULL,
                    pictures = TRUE, example = FALSE, filter = NULL, peakMonitor = NULL, pic_extension = c(".tiff",
                    ".png"), parallel = NULL, group = NULL, derivatization = NULL, cores = 1,
                    column_set = NULL, prog = NULL, ion_mode = NULL, plot_eic = NULL,
                    RI = NULL, replicate = NULL, mergeCompounds = NULL, removeCompounds = NULL,
                    info = NULL, Ri = NULL, pre_anno = NULL, min_peaks = 5, sample_names = NULL) {

  # ===============================================================
  # DEFINIÇÃO DE AMBIENTE: EXEMPLO OU USUÁRIO (UNIFICADO)
  # ===============================================================
  if (example == TRUE) {
    message("Example mode activated. Loading built-in data...")

    # 1. Busca o CSV de metadados dentro do pacote
    meta_path <- system.file("extdata", "metadata.csv", package = "PipMet")

    if (meta_path == "") {
      stop("Example metadata.csv not found in inst/extdata!")
    }

    # 2. Lê o CSV (Evita o erro de 'object metadata not found' do RData)
    metadata <- read.csv(meta_path, stringsAsFactors = FALSE)

    # 3. Ajusta caminhos e garante extensões nos nomes dos arquivos
    metadata$file <- sapply(metadata$file, function(x) {
      file_name <- basename(x)

      # Garante que tenha .mzXML (ajuste conforme os arquivos da sua pasta)
      if (!grepl("\\.mzXML$|\\.mzML$", file_name, ignore.case = TRUE)) {
        file_name <- paste0(file_name, ".mzXML")
      }

      system.file("extdata", file_name, package = "PipMet")
    })

    # 4. Configurações de diretório temporário para o exemplo
    myDir <- tempdir()
    sample_dir <- system.file("extdata", package = "PipMet")
    extension <- ".mzXML"

  } else {
    # --- MODO NORMAL (FALSE): Seleção Interativa Única ---
    message("Interactive mode: Please select your files and folders.")

    if (is.null(myDir)) {
      myDir <- svDialogs::dlg_dir(default = getwd(), title = "Select Working Directory")$res
    }

    if (is.null(sample_dir)) {
      sample_dir <- svDialogs::dlg_dir(default = myDir, title = "Select Samples Directory")$res
    }

    if (is.null(metadata)) {
      meta_res <- svDialogs::dlg_open(default = myDir, title = "Select Metadata File",
                                      filters = svDialogs::dlg_filters[c("csv"), ])$res
      if (length(meta_res) > 0) {
        metadata <- read.csv(meta_res, stringsAsFactors = FALSE)

        # Garante caminhos absolutos para o XCMS não se perder
        metadata$file <- file.path(sample_dir, basename(metadata$file))
      }
    }

    # Validação de segurança: se o usuário cancelar as janelas
    if (length(myDir) == 0 || length(sample_dir) == 0 || is.null(metadata)) {
      stop("Process cancelled by user: Paths or metadata not provided.")
    }
  }
  # ===============================================================
  # --- 1. CONFIGURAÇÕES E MONITORIZAÇÃO DE IÕES ---
  if (pictures == TRUE) {
    if (is.null(peakMonitor)) {
      peakMonitor <- svDialogs::dlg_list(c("Yes", "No"), multiple = FALSE, title = "Monitor EICs?")$res
    }
    ions <- list()
    if (peakMonitor == "Yes" | peakMonitor == TRUE) {
      okay <- 1; ei <- 1
      while (okay == 1) {
        ions[[ei]] <- vector(mode = "list", length = 2)
        names(ions[[ei]]) <- c("mz", "rt")
        ions[[ei]][["mz"]] <- as.integer(svDialogs::dlgInput(paste0("Mz of EIC ", ei, ":"), "0")$res)
        ions[[ei]][["rt"]] <- as.integer(svDialogs::dlgInput(paste0("Rt of EIC ", ei, " (+/- 5s):"), "0")$res)
        okay <- menu(c("Yes", "No"), graphics = TRUE, title = "Monitor another one?")
        ei <- ei + 1
      }
    }
  }

  # --- 2. PARALELIZAÇÃO ---
  if (is.null(parallel)) {
    parallel <- svDialogs::dlg_list(c("Serial Param", "Snow Param", "MultiCore Param"),
                                 multiple = FALSE, title = "Parallelization mode:")$res
  }
  if (parallel == "Serial Param") { register(SerialParam(), default = TRUE) }
  if (parallel == "Snow Param") {
    if (is.null(cores)) {
      cores <- svDialogs::dlgInput(paste0("Number of cores (", detectCores(), " available):"), detectCores() - 2)$res
    }
    register(SnowParam(workers = as.numeric(cores)), default = TRUE)
  }
  if (parallel == "MultiCore Param") { register(MulticoreParam(), default = TRUE) }

  # --- 3. LEITURA E PRÉ-PROCESSAMENTO ---
  message("Reading files...")
  quiet(read <- read_data(peakMonitor = peakMonitor, ions = ions, sample_dir = sample_dir,
                          metadata = metadata, extension = extension, myDir = myDir, pictures = pictures,
                          example = example, pic_extension = pic_extension, sample_names = sample_names))
  colors <- read[[1]]; metadata <- read[[2]]; raw_data <- read[[3]]; myDir <- read[[4]]

  message("Pre-processing files...")
  quiet(processed <- process(raw_data = raw_data, metadata = metadata, myDir = myDir,
                               colors = colors, peakMonitor = peakMonitor, ions = ions, pictures = pictures,
                               filter = filter, pic_extension = pic_extension, group = group))
  xdata4 <- processed$xdata4; group <- processed$group

  # --- 4. AGRUPAMENTO E ANOTAÇÃO (NIST) ---
  message("Grouping peaks into spectra...")
  quiet(spectra <- getSpectra(xdata4, raw_data, min_peaks, colors, column_set = column_set,
                               prog = prog, ion_mode = ion_mode, plot_eic = plot_eic))
  anIC <- spectra[[1]]; pslist <- spectra[[2]]; result <- spectra[[3]]

  svDialogs::dlg_message("Annotatation step: The files 'pre_anno.csv' and 'spectra.msp' were created in you directory. Upload the file 'spectra.msp' in NIST MS Search and annotate the spectra in the file 'pre_anno', in the column 'Annotation', according to the spectra 'id'. After, press 'ok'.")$res

  calculateRI(RI, result)

  message("Annotating spectra ...")
  quiet(annot <- annot_images(pslist, myDir, pictures, pic_extension = pic_extension,
                                pre_anno = pre_anno, example = example))
  apslist <- annot$apslist; pre_anno <- annot$r

  # --- 5. NORMALIZAÇÃO ---
  message("Normalizing data...")
  quiet(n <- normalize_data(anIC, pslist, metadata, myDir, pre_anno, pic_extension = pic_extension,
                            derivatization, mergeCompounds, removeCompounds, group))

  # --- 6. STATISTICS AND FIGURES (FULL INTERACTIVE CONTROL) ---
  if (nrow(metadata) > 1 & nrow(pre_anno) > 1) {
    if (is.null(replicate)) {
      if (!"tec_rep" %in% colnames(metadata)) {
        replicate <- svDialogs::dlg_list(c(colnames(metadata), "No information"),
                                      multiple = FALSE, title = "Replicate column:")$res
      }
    } else { replicate <- "tec_rep" }

    if (pictures == TRUE) {
      # Initial selection filter
      pic <- svDialogs::dlg_list(c("Volcano - level 1", "Volcano - level 2", "PCA Analysis", "Heatmaps", "All"),
                                multiple = TRUE, title = "Initial selection of Statistics pictures:")$res

      # --- VOLCANO LEVEL 1 ---
      if ("Volcano - level 1" %in% pic | "All" %in% pic) {
        if (askYesNo("Would you like to perform Volcano Level 1?") == TRUE) {
          okay_v1 <- 1
          while (okay_v1 == 1) {
            setwd(myDir); try(vol_lvl1(n, metadata, myDir, pic_extension = pic_extension))
            okay_v1 <- menu(c("Repeat Level 1", "Next Step"), graphics = TRUE, title = "Volcano Level 1 Done")
          }
        }
      }

      # --- VOLCANO LEVEL 2 ---
      if ("Volcano - level 2" %in% pic | "All" %in% pic) {
        if (askYesNo("Would you like to perform Volcano Level 2 (Interaction Analysis)?") == TRUE) {
          okay_v2 <- 1
          while (okay_v2 == 1) {
            if (!exists("volDir")) { volDir <- "Statistics/Volcanos" }
            setwd(myDir); try(vol_lvl2(n, metadata, myDir, volDir, pic_extension = pic_extension))
            okay_v2 <- menu(c("Repeat Level 2", "Next Step (PCA)"), graphics = TRUE, title = "Volcano Level 2 Done")
          }
        }
      }

      # --- PCA ANALYSIS (GROUPS & REPETITION) ---
      if ("PCA Analysis" %in% pic | "All" %in% pic) {
        if (askYesNo("Would you like to perform PCA Analysis?") == TRUE) {
          okay_pca <- 1
          while (okay_pca == 1) {

            # Sub-menu for PCA type
            pca_type <- menu(c("GENERAL PCA (All metabolites)", "IDENTIFIED COMPOUNDS PCA (Annotated Only)"),
                             graphics = TRUE, title = "Which type of PCA would you like to run now?")

            if (pca_type == 1) {
              message("Generating General PCA...")
              setwd(myDir); try(PCA_general(n, metadata, myDir, colors, pic_extension = pic_extension))
            } else if (pca_type == 2) {
              # Filter for identified compounds
              n_id <- n[!is.na(n[, 3]) & n[, 3] != "", ]
              if (nrow(n_id) > 2) {
                message("Generating PCA for Identified Compounds...")
                setwd(myDir); try(PCA_identified(n_id, metadata, myDir, colors, pic_extension = pic_extension))
              } else { svDialogs::dlg_message("Not enough identified compounds to generate PCA.")$res }
            }

            # Choice to repeat for other groups or types
            okay_pca <- menu(c("Repeat PCA (Select different groups or types)", "Next Step (Heatmaps)"),
                             graphics = TRUE, title = "PCA Finished. What's next?")
          }
        }
      }

      # --- HEATMAPS ---
      if ("Heatmaps" %in% pic | "All" %in% pic) {
        if (askYesNo("Would you like to generate Heatmaps?") == TRUE) {
          okay_hm <- 1
          while (okay_hm == 1) {
            message("Generating Heatmaps...")
            setwd(myDir); try(heatmap(n, metadata, myDir, colors, pic_extension = pic_extension, replicate = replicate))

            # Choice to repeat for other groups
            okay_hm <- menu(c("Repeat Heatmap (Different groups)", "Finish Statistics"),
                            graphics = TRUE, title = "Heatmap Done. Repeat or Finish?")
          }
        }
      }
    }
  }

  svDialogs::dlg_message("Processing done!")$res
# Retorno formal dos resultados (Essencial para o funcionamento do pacote e evitar erros)
  return(list(myDir = myDir, 
              metadata = metadata, 
              apslist = apslist, 
              raw_data = raw_data,
              xdata4 = xdata4, 
              anIC = anIC, 
              colors = colors, 
              quantification_table = n))
}