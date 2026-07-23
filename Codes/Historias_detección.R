# ==============================================================================
# HISTORIAS DE DETECCIÓN — BIOACÚSTICA FPVA
# Objetivo: visualizar la matriz de detección (plot(umf)) para todas las
#           especies con ventanas de 2 a 7 días, en mosaico por especie
# ==============================================================================

library(unmarked)
library(lubridate)
library(dplyr)
library(tidyr)
library(openxlsx)
library(data.table)
library(ggplot2)
library(patchwork)

# ==========================================================
# 1. CONFIGURACIÓN
# ==========================================================
path_data <- "~/Desktop/FPVA/Data/Bioacustica/Plantilla monitoreo acústico - FPV Amazonía.xlsx"
path_out  <- "~/Desktop/FPVA/Resultados/Historias_Deteccion/"
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

obs_data <- fread("~/Desktop/FPVA/Data/Bioacustica/observations.csv")
dev_data <- fread("~/Desktop/FPVA/Data/Bioacustica/deployments.csv")
med_data <- read.xlsx(path_data, sheet = "Media")

ventanas <- 2:7

especies_objetivo <- c(
  "Cyanocorax violaceus",
  "Myrmothera campanisona",
  "Myrmoborus myotherinus",
  "Ramphastos tucanus",
  "Lipaugus vociferans",
  "Akletos melanoceps",
  "Hylophylax naevius",
  "Thamnomanes ardesiacus",
  "Thamnomanes caesius"
)

# ==========================================================
# 2. PREPARAR DETECCIONES
# ==========================================================
obs_con_tiempo <- obs_data |>
  left_join(med_data |> select(mediaID, timestamp), by = "mediaID") |>
  mutate(time_obj = ymd_hms(timestamp)) |>
  filter(!is.na(time_obj))

fecha_min <- min(as.Date(ymd_hms(dev_data$deploymentStart)), na.rm = TRUE)
fecha_max <- max(as.Date(ymd_hms(dev_data$deploymentEnd)),   na.rm = TRUE)
ids_todos <- dev_data$deploymentID

cat("Marco temporal:", as.character(fecha_min), "->", as.character(fecha_max), "\n")
cat("Sitios:", length(ids_todos), "\n\n")

# ==========================================================
# 3. FUNCIÓN: construir y_matrix para una especie y ventana
# ==========================================================
construir_y <- function(especie, ventana) {
  
  n_oc <- as.integer(ceiling(
    as.numeric(difftime(fecha_max, fecha_min, units = "days")) / ventana))
  
  obs_sp <- obs_con_tiempo |>
    filter(scientificName == especie) |>
    mutate(fecha   = as.Date(time_obj),
           ocasion = as.integer(floor(
             as.numeric(difftime(fecha, fecha_min, units = "days")) / ventana) + 1)) |>
    filter(ocasion >= 1, ocasion <= n_oc)
  
  y_mat <- obs_sp |>
    group_by(deploymentID, ocasion) |>
    summarise(presencia = 1L, .groups = "drop") |>
    complete(deploymentID = ids_todos,
             ocasion      = 1:n_oc,
             fill         = list(presencia = 0L)) |>
    pivot_wider(names_from  = ocasion,
                values_from = presencia,
                values_fill = 0L,
                names_prefix = "oc_") |>
    arrange(match(deploymentID, ids_todos))
  
  y_final <- y_mat |> select(starts_with("oc_")) |> as.matrix()
  y_final[is.na(y_final)] <- 0L
  
  list(
    y         = y_final,
    n_oc      = n_oc,
    n_det     = sum(y_final),
    n_sitios  = sum(rowSums(y_final) > 0),
    prop_occ  = round(mean(rowSums(y_final) > 0), 3)
  )
}

# ==========================================================
# 4. FUNCIÓN: graficar historia de detección como heatmap
# ==========================================================
plot_historia <- function(especie, ventana, y_info) {
  
  y   <- y_info$y
  slug <- gsub(" ", "_", especie)
  nombre_corto <- word(especie, 2)  # solo epiteto
  
  # Convertir a data.frame largo para ggplot
  df <- as.data.frame(y)
  df$sitio <- seq_len(nrow(df))
  df_long <- df |>
    pivot_longer(-sitio, names_to = "ocasion", values_to = "deteccion") |>
    mutate(ocasion = as.integer(str_remove(ocasion, "oc_")),
           deteccion = factor(deteccion, levels = c(0, 1),
                              labels = c("Ausente", "Presente")))
  
  # Ordenar sitios por número de detecciones
  orden_sitios <- df |>
    mutate(total = rowSums(across(-sitio))) |>
    arrange(desc(total)) |>
    pull(sitio)
  df_long$sitio <- factor(df_long$sitio, levels = orden_sitios)
  
  titulo <- paste0(nombre_corto, " | v=", ventana, "d")
  subtit <- paste0("Det:", y_info$n_det,
                   " | Sit:", y_info$n_sitios,
                   "/", nrow(y),
                   " (", round(y_info$prop_occ*100), "%)")
  
  ggplot(df_long, aes(x = ocasion, y = sitio, fill = deteccion)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_manual(values = c("Ausente" = "#f0f0f0", "Presente" = "#2166ac"),
                      guide = "none") +
    scale_x_continuous(breaks = seq(1, max(df_long$ocasion), by = 2),
                       expand = c(0, 0)) +
    labs(title    = titulo,
         subtitle = subtit,
         x = "Ocasion", y = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      plot.title    = element_text(face = "bold.italic", size = 7),
      plot.subtitle = element_text(size = 6, color = "grey40"),
      axis.text.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      panel.grid    = element_blank(),
      plot.margin   = margin(3, 3, 3, 3)
    )
}

# ==========================================================
# 5. GENERAR TODOS LOS PLOTS Y MOSAICOS
# ==========================================================
cat("=== GENERANDO HISTORIAS DE DETECCIÓN ===\n\n")

# Lista para guardar todos los plots
todos_plots <- list()

for (especie in especies_objetivo) {
  cat("Especie:", especie, "\n")
  plots_especie <- list()
  
  for (v in ventanas) {
    y_info <- tryCatch(
      construir_y(especie, v),
      error = function(e) {
        cat("  ERROR ventana", v, ":", conditionMessage(e), "\n")
        NULL
      }
    )
    
    if (!is.null(y_info)) {
      cat(sprintf("  ventana=%dd | det=%d | sitios=%d/%d\n",
                  v, y_info$n_det, y_info$n_sitios, length(ids_todos)))
      
      p <- plot_historia(especie, v, y_info)
      plots_especie[[as.character(v)]] <- p
    }
  }
  
  # Mosaico por especie (6 ventanas en 2 filas x 3 columnas)
  if (length(plots_especie) > 0) {
    slug  <- gsub(" ", "_", especie)
    panel <- wrap_plots(plots_especie, ncol = 3) +
      plot_annotation(
        title    = paste("Historia de detección —", especie),
        subtitle = paste0("Ventanas de ", min(ventanas), " a ",
                          max(ventanas), " días | Azul = detección"),
        theme = theme(
          plot.title    = element_text(size = 11, face = "bold.italic"),
          plot.subtitle = element_text(size = 9, color = "grey40")
        )
      )
    
    ggsave(file.path(path_out, paste0("Historia_Deteccion_", slug, ".png")),
           plot = panel, width = 12, height = 7, dpi = 200)
    
    todos_plots[[especie]] <- plots_especie
    cat("  Panel guardado.\n")
  }
  cat("\n")
}

# ==========================================================
# 6. MOSAICO GLOBAL — todas las especies x todas las ventanas
# ==========================================================
cat("Generando mosaico global...\n")

# Seleccionar ventana 3 y 6 para el mosaico global (más compacto)
ventanas_resumen <- c(3, 6)
plots_global <- list()

for (especie in especies_objetivo) {
  for (v in ventanas_resumen) {
    y_info <- tryCatch(construir_y(especie, v), error = function(e) NULL)
    if (!is.null(y_info)) {
      p <- plot_historia(especie, v, y_info)
      plots_global[[paste(especie, v, sep="|")]] <- p
    }
  }
}

if (length(plots_global) > 0) {
  # 9 especies x 2 ventanas = 18 paneles, 6 columnas
  mosaico_global <- wrap_plots(plots_global, ncol = 6) +
    plot_annotation(
      title    = "Historias de detección — Bioacústica FPVA",
      subtitle = "Ventanas de 3 y 6 días | Azul = detección | Ordenado por n° detecciones",
      caption  = paste0("Marco temporal: ", fecha_min, " → ", fecha_max,
                        " | ", length(ids_todos), " sitios"),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        plot.caption  = element_text(size = 8,  color = "grey60")
      )
    )
  
  ggsave(file.path(path_out, "Mosaico_Global_Bioacustica.png"),
         plot = mosaico_global, width = 20, height = 12, dpi = 200)
  cat("Mosaico global guardado.\n")
}

# ==========================================================
# 7. TABLA RESUMEN DE DETECCIONES
# ==========================================================
cat("\nGenerando tabla resumen...\n")

tabla_det <- expand.grid(
  Especie      = especies_objetivo,
  Ventana_dias = ventanas,
  stringsAsFactors = FALSE
) |>
  rowwise() |>
  mutate(
    info = list(tryCatch(construir_y(Especie, Ventana_dias), error=function(e) NULL)),
    N_ocasiones    = if (!is.null(info)) info$n_oc     else NA,
    N_detecciones  = if (!is.null(info)) info$n_det    else NA,
    N_sitios_detec = if (!is.null(info)) info$n_sitios else NA,
    Prop_sitios    = if (!is.null(info)) info$prop_occ else NA
  ) |>
  ungroup() |>
  select(-info) |>
  mutate(Especie_corta = word(Especie, 2)) |>
  arrange(Especie, Ventana_dias)

write.csv(tabla_det,
          file.path(path_out, "Tabla_Detecciones_por_Ventana.csv"),
          row.names = FALSE)

cat("\n=== TABLA RESUMEN ===\n")
print(tabla_det |> select(Especie_corta, Ventana_dias,
                          N_ocasiones, N_detecciones,
                          N_sitios_detec, Prop_sitios),
      n = Inf)

cat("\nResultados guardados en:", path_out, "\n")
cat("=== COMPLETADO ===\n")