#' Statistics figures - Volcano Level 2 (Flexible Interaction)
#'
#' This function plots volcanos with two levels of comparison, allowing independent
#' selection of target and baseline for two different pairs.
#' @export
vol_lvl2 <- function(n, metadata, myDir, volDir, pic_extension = c(".tiff", ".png")) {
  
  # 1. Preparação da Matriz de Dados
  mat <- n[, 6:ncol(n)]
  mat <- apply(mat, MARGIN = 2, FUN = as.numeric)
  rownames(mat) <- n[, 1]
  mat[is.na(mat)] <- 0

  # 2. Interface de Seleção Flexível
  dlg_message("Volcano Level 2: Define two pairs (Target vs Baseline) to compare their responses.")$res

  col_idx <- menu(names(metadata), graphics = TRUE, title = "Select the Metadata Column (Group/Condition)")
  if (col_idx == 0) return(message("Cancelled."))
  
  all_lvls <- as.character(unique(metadata[, col_idx]))

  # Seleção do Par 1
  dlg_message("PAIR 1: Select the TARGET then the BASELINE.")$res
  p1_t <- select.list(all_lvls, title = "Pair 1: Target", graphics = TRUE)
  p1_b <- select.list(all_lvls, title = "Pair 1: Baseline", graphics = TRUE)
  
  # Seleção do Par 2
  dlg_message("PAIR 2: Select the TARGET then the BASELINE.")$res
  p2_t <- select.list(all_lvls, title = "Pair 2: Target", graphics = TRUE)
  p2_b <- select.list(all_lvls, title = "Pair 2: Baseline", graphics = TRUE)

  if (any(c(p1_t, p1_b, p2_t, p2_b) == "")) return(message("Selection incomplete."))

  # 3. CÁLCULO ESTATÍSTICO DE INTERAÇÃO
  idx1t <- which(metadata[, col_idx] == p1_t)
  idx1b <- which(metadata[, col_idx] == p1_b)
  idx2t <- which(metadata[, col_idx] == p2_t)
  idx2b <- which(metadata[, col_idx] == p2_b)

  interaction_results <- sapply(1:nrow(mat), function(i) {
    resp1 <- mat[i, idx1t] - mean(mat[i, idx1b])
    resp2 <- mat[i, idx2t] - mean(mat[i, idx2b])
    delta_fc <- mean(resp1) - mean(resp2)
    p_val <- tryCatch({
      t.test(resp1, resp2)$p.value
    }, error = function(e) return(1.0))
    return(c(delta_fc = delta_fc, p_val = p_val))
  })

  # 4. Criação da Tabela Reporte (Mantendo colunas da imagem)
  de <- data.frame(
    id = n[, 1],                         # ID original (crucial para rastreio)
    Fragment_Ion = n[, 2],               # Informação técnica
    Compound_Name = ifelse(is.na(n[, 3]) | n[, 3] == "", n[, 1], n[, 3]),
    Metabolic_Class = n[, 5],            # Classificação química
    Delta_Log2FC = interaction_results["delta_fc", ],
    P_Interaction = interaction_results["p_val", ]
  )

  # Classificação para o Plot e Tabela
  label_p1 <- paste0("Higher response in ", p1_t, "/", p1_b)
  label_p2 <- paste0("Higher response in ", p2_t, "/", p2_b)
  
  de$Status <- "Consistent Response"
  de$Status[de$Delta_Log2FC > 0.6 & de$P_Interaction < 0.05] <- label_p1
  de$Status[de$Delta_Log2FC < -0.6 & de$P_Interaction < 0.05] <- label_p2
  de$delabel <- ifelse(de$Status != "Consistent Response", de$Compound_Name, NA)

  # 5. Geração do Gráfico
  p <- ggplot(data = de, aes(x = Delta_Log2FC, y = -log10(P_Interaction), col = Status, label = delabel)) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual(values = setNames(c("red", "gray70", "blue"), c(label_p1, "Consistent Response", label_p2))) +
    theme_minimal() +
    geom_vline(xintercept = c(-0.6, 0.6), col = "black", linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", alpha = 0.3) +
    labs(title = "Level 2: Differential Response Analysis",
         subtitle = paste0("Contrast: (", p1_t, "-", p1_b, ") vs (", p2_t, "-", p2_b, ")"),
         x = "Delta Fold Change (Difference of Responses)",
         y = "-Log10 p-value (Interaction Significance)") +
    geom_text_repel(max.overlaps = 20) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

  # 6. Salvamento Organizado
  main_path <- file.path(volDir, "Volcano_Level2_Final")
  if (!dir.exists(main_path)) dir.create(main_path, recursive = TRUE)
  
  for (ext in pic_extension) {
    ggsave(file.path(main_path, paste0("Volcano_Level2", ext)), plot = p, width = 18, height = 16, units = "cm", dpi = 300)
  }

  write.csv(de, file = file.path(main_path, "Significance_Level2_Table.csv"), row.names = FALSE)
  
  setwd(myDir)
  dlg_message("Level 2 analysis complete. The table includes IDs, Names and Metabolic Classes.")$res
}