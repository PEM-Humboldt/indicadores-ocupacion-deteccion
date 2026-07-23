# ==============================================================================
# ANÁLISIS DE OCUPACIÓN - FOTOTRAMPEO - FPVA
# Especies objetivo: Tapirus terrestris, Cuniculus paca, Dasypus novemcinctus,
#                    Dicotyles tajacu, Tayassu pecari, Odocoileus virginianus,
#                    Nasua nasua
#
# Metodología: Modelo de ocupación de MacKenzie et al. (2002) via `unmarked`
# Selección de modelos: dredge() basado en buenas prácticas eBird/GBIF
# Referencia: https://strimas.com/ebp-workshop-au/occupancy.html
# ==============================================================================

# ==========================================================
# 0. LIBRERÍAS
# ==========================================================
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(unmarked)
library(MuMIn)       # Para dredge() y model.avg()
library(ggplot2)
library(lubridate)
library(tidyterra)
library(openxlsx)
library(patchwork)

# ==========================================================
# 1. CONFIGURACIÓN GENERAL
# ==========================================================

path_data <- "~/Desktop/FPVA/Data/Fototrampeo/I2D_FPVA_Fototrampeo_20260219.xlsx"
path_out  <- "~/Desktop/FPVA/Resultados/Fototrampeo/"
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

obs_data <- read.xlsx(path_data, sheet = "Observations")
med_data <- read.xlsx(path_data, sheet = "Media")
dev_data <- read.xlsx(path_data, sheet = "Deployment")

# Vector de especies objetivo
especies_objetivo <- c(
  "Tapirus terrestris",
  "Cuniculus paca",
  "Dasypus novemcinctus",
  "Dicotyles tajacu",
  "Tayassu pecari",
  "Odocoileus virginianus",
  "Nasua nasua"
)

# ==========================================================
# 2. PREPARACIÓN DE COVARIABLES (Limpieza y Estandarización)
# ==========================================================
# --- 2.1 Covariables originales ---
path_base_orig <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/Co-variables Ocupación/"

cov_stack_orig <- c(
  rast(paste0(path_base_orig, "BosqueNoBosque__90mts/agreg_bosque90_2.tif")),
  rast(paste0(path_base_orig, "EucDistance_Drenajes90mts/RiosEudis90.tif")),
  rast(paste0(path_base_orig, "EucliDistance_vias90mts/ViasEudis90.tif"))
)
names(cov_stack_orig) <- c("bosque", "dist_rios", "dist_vias")

# --- 2.2 Nuevas covariables ---
path_base_new <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/VariablesOcupacion_2/"

r_ndvi  <- rast(paste0(path_base_new, "indicesNDVI_WI/ndvi_300725.tif"))
r_ndwi  <- rast(paste0(path_base_new, "indicesNDVI_WI/ndwi_300725.tif"))
r_d_com <- rast(paste0(path_base_new, "DrenajesCompleto_2.tif/DrenajesCompleto_2.tif"))
r_d_pri <- rast(paste0(path_base_new, "DrenajesPrincip_2.tif/DrenajesPrincip_2.tif"))
r_d_sec <- rast(paste0(path_base_new, "DrenajesSecund_2.tif/DrenajesSecund_2.tif"))
r_v_com <- rast(paste0(path_base_new, "viascompletas_2/viascompletas_2.tif"))
r_v_sec <- rast(paste0(path_base_new, "viassecundarias_2/viassecundarias_2.tif"))
r_v_pri <- rast(paste0(path_base_new, "Viasprincipales_2.tif/Viasprincipales_2.tif"))
r_casas <- rast(paste0(path_base_new, "equipamien_eudist_2.tif/equipamien_eudist_2.tif"))

# ==========================================================
# 2.3 Alinear todas las capas al mismo grid de referencia
#     antes de construir el stack
# ==========================================================

# La capa de bosque será la referencia (define extent, resolución y CRS)
r_ref <- cov_stack_orig[[1]]  # bosque

# Función auxiliar: reproyectar + recortar + remuestrear al grid de referencia
alinear_capa <- function(r, referencia) {
  # 1. Si el CRS es distinto, reproyectar primero
  if (!same.crs(r, referencia)) {
    r <- project(r, crs(referencia))
  }
  # 2. Recortar al extent de referencia (con tolerancia)
  r <- crop(r, ext(referencia))
  # 3. Remuestrear al mismo grid exacto (resolución + origen)
  r <- resample(r, referencia, method = "bilinear")
  return(r)
}

# Alinear covariables originales entre sí
r_bosque   <- cov_stack_orig[[1]]  # ya es la referencia
r_dist_rio <- alinear_capa(cov_stack_orig[[2]], r_ref)
r_dist_via <- alinear_capa(cov_stack_orig[[3]], r_ref)

# Alinear nuevas covariables
r_ndvi_a     <- alinear_capa(r_ndvi,  r_ref)
r_ndwi_a     <- alinear_capa(r_ndwi,  r_ref)
r_d_com_a    <- alinear_capa(r_d_com, r_ref)
r_d_pri_a    <- alinear_capa(r_d_pri, r_ref)
r_d_sec_a    <- alinear_capa(r_d_sec, r_ref)
r_v_com_a    <- alinear_capa(r_v_com, r_ref)
r_v_sec_a    <- alinear_capa(r_v_sec, r_ref)
r_v_pri_a    <- alinear_capa(r_v_pri, r_ref)
r_casas_a    <- alinear_capa(r_casas, r_ref)

# Construir el stack ya alineado
cov_stack <- c(
  r_bosque, r_dist_rio, r_dist_via,
  r_ndvi_a, r_ndwi_a,
  r_d_com_a, r_d_pri_a, r_d_sec_a,
  r_v_com_a, r_v_sec_a, r_v_pri_a,
  r_casas_a
)
names(cov_stack) <- c(
  "bosque", "dist_rios", "dist_vias",
  "ndvi", "ndwi",
  "dist_dren_com", "dist_dren_pri", "dist_dren_sec",
  "dist_vias_com", "dist_vias_sec", "dist_vias_pri",
  "dist_casas"
)

# Verificación rápida — debe imprimir TRUE
cat("¿Stack alineado correctamente?", compareGeom(cov_stack[[1]], cov_stack[[12]], 
                                                  stopOnError = FALSE), "\n")

# --- 2.4 Reproyección de todas las capas a la misma referencia ---
# Reproyectamos al CRS de la primera capa (bosque) como referencia
ref_crs <- crs(cov_stack$bosque)
cov_stack <- lapply(names(cov_stack), function(nm) {
  r <- cov_stack[[nm]]
  if (!compareGeom(r, cov_stack$bosque, stopOnError = FALSE)) {
    r <- project(r, ref_crs)
    r <- resample(r, cov_stack$bosque, method = "bilinear")
  }
  r
}) |> do.call(what = c)
names(cov_stack) <- c(
  "bosque", "dist_rios", "dist_vias",
  "ndvi", "ndwi",
  "dist_dren_com", "dist_dren_pri", "dist_dren_sec",
  "dist_vias_com", "dist_vias_sec", "dist_vias_pri",
  "dist_casas"
)

# --- 2.5 Extracción en puntos de fototrampeo ---
puntos_sf <- dev_data |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(cov_stack))

site_covs_raw <- terra::extract(cov_stack, puntos_sf)
site_covs_raw$deploymentID <- dev_data$deploymentID

# Nombres de todas las covariables de sitio (sin ID)
cov_names <- names(cov_stack)

# Guardar parámetros de escalado (media y SD) — necesarios para mapas
params <- site_covs_raw |>
  summarise(across(all_of(cov_names),
                   list(mu = ~mean(.x, na.rm = TRUE),
                        sd  = ~sd(.x,   na.rm = TRUE))))

site_covs_scaled <- site_covs_raw |>
  mutate(across(all_of(cov_names), scale)) |>
  filter(complete.cases(pick(all_of(cov_names))))

# ==========================================================
# 3. FUNCIÓN PRINCIPAL DE ANÁLISIS POR ESPECIE
# ==========================================================
analizar_especie <- function(especie_objetivo) {
  
  cat("\n========================================================\n")
  cat("  Procesando:", especie_objetivo, "\n")
  cat("========================================================\n")
  
  # Slug para nombres de archivos
  slug <- gsub(" ", "_", especie_objetivo)
  
  # Carpeta de salida por especie
  path_spp <- file.path(path_out, slug)
  if (!dir.exists(path_spp)) dir.create(path_spp, recursive = TRUE)
  
  # --------------------------------------------------------
  # 3.1 Historia de detección
  # --------------------------------------------------------
  
  obs_con_tiempo <- obs_data |>
    left_join(med_data |> select(mediaID, timestamp), by = "mediaID") |>
    mutate(time_obj = ymd_hms(timestamp)) |>
    filter(!is.na(time_obj))
  
  # CORRECCIÓN 1: fecha_min y fecha_max desde el despliegue de cámaras,
  # no desde las observaciones. Así todos los sitios comparten el mismo
  # marco temporal independientemente de cuándo detectaron cada especie.
  fecha_min <- min(as.Date(ymd_hms(dev_data$deploymentStart)), na.rm = TRUE)
  fecha_max <- max(as.Date(ymd_hms(dev_data$deploymentEnd)),   na.rm = TRUE)
  
  ventana_dias <- 7   # tamaño de cada ocasión en días
  n_ocasiones_total <- as.integer(
    ceiling(as.numeric(difftime(fecha_max, fecha_min, units = "days")) / ventana_dias)
  )
  
  cat("  Ventana temporal del estudio:", as.character(fecha_min),
      "→", as.character(fecha_max), "\n")
  cat("  Número total de ocasiones (", ventana_dias, "días c/u):",
      n_ocasiones_total, "\n")
  
  # CORRECCIÓN 2: filtrar observaciones de la especie objetivo
  # y asignar ocasión dentro del marco temporal del estudio
  obs_sp <- obs_con_tiempo |>
    filter(scientificName == especie_objetivo) |>
    mutate(
      fecha   = as.Date(time_obj),
      ocasion = as.integer(
        floor(as.numeric(difftime(fecha, fecha_min, units = "days")) / ventana_dias) + 1
      )
    ) |>
    # Descartar registros fuera del marco temporal (antes o después del estudio)
    filter(ocasion >= 1, ocasion <= n_ocasiones_total)
  
  cat("  Registros de", especie_objetivo, "dentro del marco temporal:",
      nrow(obs_sp), "\n")
  
  # CORRECCIÓN 3: construir matriz con todas las ocasiones como columnas,
  # incluso las que no tienen detección (valores fijos)
  todas_ocasiones <- 1:n_ocasiones_total
  
  y_matrix <- obs_sp |>
    group_by(deploymentID, ocasion) |>
    summarise(presencia = 1L, .groups = "drop") |>
    # Completar todas las combinaciones sitio x ocasión
    complete(
      deploymentID = site_covs_scaled$deploymentID,
      ocasion      = todas_ocasiones,
      fill         = list(presencia = 0L)
    ) |>
    pivot_wider(names_from  = ocasion,
                values_from = presencia,
                values_fill = 0L,
                names_prefix = "oc_")   # prefijo para evitar columnas numéricas sueltas
  
  # CORRECCIÓN 4: alinear al orden de sitios en site_covs_scaled
  y_final <- site_covs_scaled |>
    select(deploymentID) |>
    left_join(y_matrix, by = "deploymentID") |>
    select(starts_with("oc_")) |>
    as.matrix()
  y_final[is.na(y_final)] <- 0L
  
  # Diagnóstico rápido
  cat("  Sitios:", nrow(y_final),
      "| Ocasiones:", ncol(y_final),
      "| Detecciones totales:", sum(y_final), "\n")
  cat("  Sitios con ≥1 detección:", sum(rowSums(y_final) > 0), "\n")
  
  if (sum(y_final, na.rm = TRUE) < 5) {
    cat("  ADVERTENCIA: Muy pocas detecciones para", especie_objetivo,
        "— análisis omitido.\n")
    return(invisible(NULL))
  }
  
  # --------------------------------------------------------
  # 3.2 Covariable de detección: hora de actividad
  # --------------------------------------------------------
  hora_sitio <- obs_con_tiempo |>
    filter(scientificName == especie_objetivo) |>
    mutate(hora = hour(time_obj)) |>
    group_by(deploymentID) |>
    summarise(hora_promedio = mean(hora, na.rm = TRUE)) |>
    right_join(data.frame(deploymentID = site_covs_scaled$deploymentID),
               by = "deploymentID") |>
    mutate(hora_promedio = replace_na(hora_promedio, mean(hora_promedio, na.rm = TRUE)))
  
  params$hora_actividad_mu <- mean(hora_sitio$hora_promedio, na.rm = TRUE)
  params$hora_actividad_sd <- sd(hora_sitio$hora_promedio,   na.rm = TRUE)
  
  # --------------------------------------------------------
  # 3.3 Esfuerzo de muestreo (matriz de obs-covs)
  # --------------------------------------------------------
  n_ocasiones <- ncol(y_final)
  effort_mat  <- matrix(15L, nrow = nrow(y_final), ncol = n_ocasiones)
  
  # --------------------------------------------------------
  # 3.4 Construcción del unmarkedFrame
  # --------------------------------------------------------
  site_covs_umf <- site_covs_scaled |>
    mutate(hora_actividad = scale(hora_sitio$hora_promedio)[, 1]) |>
    select(all_of(cov_names), hora_actividad)
  
  umf <- unmarkedFrameOccu(
    y        = y_final,
    siteCovs = site_covs_umf,
    obsCovs  = list(esfuerzo = effort_mat)
  )
  
  # --------------------------------------------------------
  # 3.5 Modelo global — con filtro automático de correlación
  # --------------------------------------------------------
  
  # CORRECCIÓN: cov_seleccionadas siempre parte siendo todas las covariables.
  # El bloque de filtro puede reducirla, pero si no hay problemas queda igual.
  cov_seleccionadas <- cov_names
  
  umbral_cor <- 0.7
  
  cor_mat <- cor(site_covs_umf[, cov_names], use = "complete.obs")
  
  repeat {
    cor_sub <- cor(site_covs_umf[, cov_seleccionadas], use = "complete.obs")
    cor_sub[lower.tri(cor_sub, diag = TRUE)] <- NA
    problema <- which(abs(cor_sub) > umbral_cor, arr.ind = TRUE)
    
    if (nrow(problema) == 0) break
    
    vars_problema <- unique(c(rownames(cor_sub)[problema[, 1]],
                              colnames(cor_sub)[problema[, 2]]))
    n_cors <- sapply(vars_problema, function(v) {
      sum(abs(cor_sub[v, ]) > umbral_cor, na.rm = TRUE) +
        sum(abs(cor_sub[, v]) > umbral_cor, na.rm = TRUE)
    })
    
    eliminar <- names(which.max(n_cors))
    cat("  Eliminando por correlación alta (|r| >", umbral_cor, "):", eliminar, "\n")
    cov_seleccionadas <- setdiff(cov_seleccionadas, eliminar)
    
    if (length(cov_seleccionadas) < 2) break
  }
  
  # Filtro por ratio eventos/parámetros
  n_det    <- sum(y_final, na.rm = TRUE)
  n_params <- length(cov_seleccionadas) + 3
  ratio_EP <- n_det / n_params
  
  cat("  Detecciones:", n_det, "| Parámetros:", n_params,
      "| Ratio E/P:", round(ratio_EP, 1), "\n")
  cat("  Covariables retenidas:", paste(cov_seleccionadas, collapse = ", "), "\n")
  
  if (ratio_EP < 10) {
    n_max_cov <- max(1, floor(n_det / 10) - 2)
    cat("  ADVERTENCIA: ratio E/P bajo. Reduciendo a", n_max_cov, "covariable(s).\n")
    prioridad <- c("bosque", "ndvi", "ndwi", "dist_rios", "dist_dren_com",
                   "dist_dren_pri", "dist_dren_sec", "dist_vias", "dist_casas",
                   "dist_vias_com", "dist_vias_sec", "dist_vias_pri")
    cov_seleccionadas <- intersect(prioridad, cov_seleccionadas)[1:n_max_cov]
    cat("  Covariables finales:", paste(cov_seleccionadas, collapse = ", "), "\n")
  }
  
  # Construir y ajustar modelo global
  formula_global <- as.formula(
    paste("~ hora_actividad + esfuerzo ~", paste(cov_seleccionadas, collapse = " + "))
  )
  cat("  Fórmula global:", deparse(formula_global), "\n")
  
  occ_global <- tryCatch(
    occu(formula_global, data = umf),
    error = function(e) {
      cat("  ERROR en modelo global:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  if (is.null(occ_global)) return(invisible(NULL))
  
  if (occ_global@opt$convergence != 0) {
    cat("  ADVERTENCIA: modelo global no convergió (código:",
        occ_global@opt$convergence, ").\n")
  }
  
  # --------------------------------------------------------
  # 3.6 SELECCIÓN DE MODELOS
  # --------------------------------------------------------
  cat("  Ejecutando dredge() — puede tardar varios minutos...\n")
  
  # Identificar términos de detección para fijarlos en todos los submodelos
  det_terms <- getAllTerms(occ_global) |>
    discard(str_detect, pattern = "psi")
  
  # Evaluamos todos los submodelos de ocupación manteniendo fija la detección
  occ_dredge <- dredge(occ_global, fixed = det_terms)
  
  # Tabla de selección
  tabla_dredge <- as.data.frame(occ_dredge) |>
    select(starts_with("psi("), df, AICc, delta, weight) |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  write.csv(tabla_dredge,
            file.path(path_spp, paste0("Tabla_Dredge_", slug, ".csv")),
            row.names = FALSE)
  
  cat("  Top 5 modelos por AICc:\n")
  print(head(tabla_dredge, 5))
  
  # --------------------------------------------------------
  # 3.7 Promediado de modelos (model averaging — top 95% peso)
  # --------------------------------------------------------
  occ_dredge_95 <- get.models(occ_dredge, subset = cumsum(weight) <= 0.95)
  occ_avg       <- model.avg(occ_dredge_95, fit = TRUE)
  
  cat("  Modelos incluidos en promediado:", length(occ_dredge_95), "\n")
  
  # Guardar resumen del modelo promediado
  sink(file.path(path_spp, paste0("Resumen_ModeloPromediado_", slug, ".txt")))
  print(summary(occ_avg))
  sink()
  
  # --------------------------------------------------------
  # 3.8 Modelo nulo para indicadores base
  # --------------------------------------------------------
  fm0 <- occu(~1 ~1, data = umf)
  
  
  # --------------------------------------------------------
  # 3.9 GRÁFICAS DE RESPUESTA — todas las covariables
  # --------------------------------------------------------
  
  # Función robusta: detecta automáticamente los nombres de columna
  # según si el modelo es model.avg o un occu simple
  get_pred_df <- function(modelo, nd, type = "state") {
    pred <- predict(modelo, newdata = nd, type = type)
    
    # model.avg devuelve: fit, se.fit (sin intervalos directos)
    # occu/predict devuelve: Predicted, SE, lower, upper
    if ("Predicted" %in% names(pred)) {
      data.frame(fit = pred$Predicted,
                 lo  = pred$lower,
                 hi  = pred$upper)
    } else {
      # Calcular IC 95% manualmente desde SE (escala logit → prob)
      data.frame(fit = pred$fit,
                 lo  = pmax(0, pred$fit - 1.96 * pred$se.fit),
                 hi  = pmin(1, pred$fit + 1.96 * pred$se.fit))
    }
  }
  
  # Función genérica de curva de respuesta
  plot_occu <- function(var_name, label_x, color_fill,
                        modelo       = occ_avg,
                        covs_en_nd   = cov_seleccionadas,
                        site_raw     = site_covs_raw) {
    
    mu_v    <- params[[paste0(var_name, "_mu")]]
    sd_v    <- params[[paste0(var_name, "_sd")]]
    
    raw_seq    <- seq(min(site_raw[[var_name]], na.rm = TRUE),
                      max(site_raw[[var_name]], na.rm = TRUE),
                      length = 200)
    scaled_seq <- (raw_seq - mu_v) / sd_v
    
    # newdata con todas las covariables del modelo en 0,
    # variando solo la variable de interés
    nd <- as.data.frame(matrix(
      0, nrow = 200,
      ncol = length(covs_en_nd) + 1,
      dimnames = list(NULL, c(covs_en_nd, "hora_actividad"))
    ))
    nd[[var_name]] <- scaled_seq
    
    pr <- get_pred_df(modelo, nd, type = "state")
    
    x_vals  <- if (var_name == "bosque") raw_seq * 100 else raw_seq
    x_label <- if (var_name == "bosque") paste(label_x, "(%)") else label_x
    
    ggplot(data.frame(x = x_vals, y = pr$fit, lo = pr$lo, hi = pr$hi),
           aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = color_fill) +
      geom_line(color = color_fill, linewidth = 1.3) +
      theme_minimal(base_size = 12) +
      labs(title    = especie_objetivo,
           subtitle = paste("Efecto de:", x_label),
           x        = x_label,
           y        = "Probabilidad de Ocupación (ψ)") +
      ylim(0, 1)
  }
  
  # Paleta y etiquetas — cubre todas las covariables posibles
  colores_occu <- c(
    bosque        = "darkgreen",
    ndvi          = "olivedrab",
    ndwi          = "steelblue",
    ndvi_ndwi     = "mediumseagreen",
    dist_rios     = "royalblue",
    dist_vias     = "firebrick",
    dist_dren_com = "dodgerblue",
    dist_dren_pri = "navy",
    dist_dren_sec = "skyblue3",
    dist_vias_com = "tomato",
    dist_vias_sec = "orangered",
    dist_vias_pri = "red4",
    dist_casas    = "sienna"
  )
  labels_occu <- c(
    bosque        = "Cobertura de Bosque",
    ndvi          = "NDVI",
    ndwi          = "NDWI",
    ndvi_ndwi     = "Índice combinado NDVI + NDWI",
    dist_rios     = "Distancia a Ríos (m)",
    dist_vias     = "Distancia a Vías (m)",
    dist_dren_com = "Distancia a Drenajes Completos (m)",
    dist_dren_pri = "Distancia a Drenajes Principales (m)",
    dist_dren_sec = "Distancia a Drenajes Secundarios (m)",
    dist_vias_com = "Distancia a Vías Completas (m)",
    dist_vias_sec = "Distancia a Vías Secundarias (m)",
    dist_vias_pri = "Distancia a Vías Principales (m)",
    dist_casas    = "Distancia a Equipamientos (m)"
  )
  
  plots_occu <- list()
  
  # --- Curva por cada covariable retenida en el modelo ---
  for (v in cov_seleccionadas) {
    cat("  Graficando:", v, "\n")
    p <- tryCatch(
      plot_occu(v, labels_occu[[v]], colores_occu[[v]]),
      error = function(e) {
        cat("  No se pudo graficar:", v, "—", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(p)) {
      plots_occu[[v]] <- p
      ggsave(
        filename = file.path(path_spp, paste0("Curva_Occu_", v, "_", slug, ".png")),
        plot = p, width = 7, height = 5, dpi = 200
      )
    }
  }
  
  # --------------------------------------------------------
  # 3.9b ÍNDICE COMBINADO NDVI + NDWI
  # --------------------------------------------------------
  if (all(c("ndvi", "ndwi") %in% cov_seleccionadas)) {
    
    cat("  Generando índice combinado NDVI + NDWI...\n")
    
    # Construir el índice: suma de Z-scores de cada índice, re-estandarizado
    ndvi_z <- as.numeric(scale(site_covs_umf$ndvi))
    ndwi_z <- as.numeric(scale(site_covs_umf$ndwi))
    indice_raw <- ndvi_z + ndwi_z
    
    site_covs_umf$ndvi_ndwi <- as.numeric(scale(indice_raw))
    
    # Parámetros para reconstruir el eje X en escala interpretable
    params$ndvi_ndwi_mu <- 0  # ya está centrado
    params$ndvi_ndwi_sd <- 1
    
    # umf y modelo con índice combinado (reemplaza ndvi y ndwi)
    cov_comb     <- c(setdiff(cov_seleccionadas, c("ndvi", "ndwi")), "ndvi_ndwi")
    formula_comb <- as.formula(
      paste("~ hora_actividad + esfuerzo ~", paste(cov_comb, collapse = " + "))
    )
    
    umf_comb <- unmarkedFrameOccu(
      y        = y_final,
      siteCovs = site_covs_umf,
      obsCovs  = list(esfuerzo = effort_mat)
    )
    
    cat("  Ajustando modelo con índice combinado:", deparse(formula_comb), "\n")
    occ_comb <- tryCatch(
      occu(formula_comb, data = umf_comb),
      error = function(e) {
        cat("  ERROR modelo combinado:", conditionMessage(e), "\n"); NULL
      }
    )
    
    if (!is.null(occ_comb)) {
      
      # Curva de respuesta del índice combinado
      idx_seq <- seq(min(site_covs_umf$ndvi_ndwi, na.rm = TRUE),
                     max(site_covs_umf$ndvi_ndwi, na.rm = TRUE),
                     length = 200)
      
      nd_comb <- as.data.frame(matrix(
        0, nrow = 200,
        ncol = length(cov_comb) + 1,
        dimnames = list(NULL, c(cov_comb, "hora_actividad"))
      ))
      nd_comb$ndvi_ndwi <- idx_seq
      
      pr_comb <- get_pred_df(occ_comb, nd_comb, type = "state")
      
      p_comb <- ggplot(
        data.frame(x = idx_seq, y = pr_comb$fit, lo = pr_comb$lo, hi = pr_comb$hi),
        aes(x = x, y = y)
      ) +
        geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25,
                    fill = colores_occu[["ndvi_ndwi"]]) +
        geom_line(color = colores_occu[["ndvi_ndwi"]], linewidth = 1.3) +
        theme_minimal(base_size = 12) +
        labs(title    = especie_objetivo,
             subtitle = "Índice combinado NDVI + NDWI",
             x        = "NDVI + NDWI (Z-scores sumados, estandarizados)",
             y        = "Probabilidad de Ocupación (ψ)",
             caption  = "Valores altos = mayor verdor y contenido hídrico") +
        ylim(0, 1)
      
      plots_occu[["ndvi_ndwi"]] <- p_comb
      ggsave(
        filename = file.path(path_spp, paste0("Curva_Occu_NDVI_NDWI_", slug, ".png")),
        plot = p_comb, width = 7, height = 5, dpi = 200
      )
      
      # Comparación AIC: combinado vs. solo NDVI vs. solo NDWI
      fm_ndvi_solo <- tryCatch(
        occu(as.formula(paste("~ hora_actividad + esfuerzo ~",
                              paste(c(setdiff(cov_seleccionadas, c("ndvi","ndwi")), "ndvi"),
                                    collapse = " + "))),
             data = umf),
        error = function(e) NULL
      )
      fm_ndwi_solo <- tryCatch(
        occu(as.formula(paste("~ hora_actividad + esfuerzo ~",
                              paste(c(setdiff(cov_seleccionadas, c("ndvi","ndwi")), "ndwi"),
                                    collapse = " + "))),
             data = umf),
        error = function(e) NULL
      )
      
      lista_comp <- list("NDVI_NDWI_combinado" = occ_comb)
      if (!is.null(fm_ndvi_solo)) lista_comp[["Solo_NDVI"]] <- fm_ndvi_solo
      if (!is.null(fm_ndwi_solo)) lista_comp[["Solo_NDWI"]] <- fm_ndwi_solo
      
      if (length(lista_comp) > 1) {
        tab_comp <- as(modSel(do.call(fitList, lista_comp)), "data.frame") |>
          mutate(across(where(is.numeric), ~ round(.x, 3)))
        cat("\n  Comparación NDVI / NDWI / combinado:\n")
        print(tab_comp)
        write.csv(tab_comp,
                  file.path(path_spp, paste0("Comparacion_NDVI_NDWI_", slug, ".csv")),
                  row.names = FALSE)
      }
    }
  } else {
    cat("  NDVI o NDWI no están en cov_seleccionadas — índice combinado omitido.\n")
  }
  
  # --------------------------------------------------------
  # 3.9c PANEL COMPLETO — todas las curvas juntas
  # --------------------------------------------------------
  plots_validos <- Filter(Negate(is.null), plots_occu)
  cat("  Total de gráficas generadas:", length(plots_validos), "\n")
  
  if (length(plots_validos) >= 2) {
    n_col  <- 3
    n_fila <- ceiling(length(plots_validos) / n_col)
    
    panel_final <- wrap_plots(plots_validos, ncol = n_col) +
      plot_annotation(
        title    = paste("Respuestas de ocupación —", especie_objetivo),
        subtitle = paste("Modelo promediado (top 95% AICc) |",
                         length(plots_validos), "covariables"),
        theme    = theme(
          plot.title    = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11, color = "grey40")
        )
      )
    
    ggsave(
      filename = file.path(path_spp, paste0("Panel_Completo_Ocupacion_", slug, ".png")),
      plot     = panel_final,
      width    = n_col * 7,
      height   = n_fila * 5,
      dpi      = 200
    )
    cat("  Panel guardado con", length(plots_validos), "gráficas.\n")
  }
  # --------------------------------------------------------
  # 3.10 MAPA DE OCUPACIÓN (con modelo promediado)
  # --------------------------------------------------------
  # Escalar el stack completo con los parámetros de los sitios de muestreo
  env_scaled <- cov_stack
  for (v in cov_names) {
    mu_v <- params[[paste0(v, "_mu")]]
    sd_v <- params[[paste0(v, "_sd")]]
    env_scaled[[v]] <- (cov_stack[[v]] - mu_v) / sd_v
  }
  
  # Superficie de predicción (todas las celdas del raster)
  pred_surface <- as.data.frame(env_scaled, xy = TRUE, na.rm = FALSE)
  pred_surface$hora_actividad <- 0  # Fijamos hora en su media (= 0 escalado)
  
  # Predicción con modelo promediado
  cat("  Generando predicciones espaciales (modelo promediado)...\n")
  occ_pred <- tryCatch(
    predict(occ_avg,
            newdata = pred_surface[, c(cov_names, "hora_actividad")],
            type    = "state"),
    error = function(e) {
      cat("  ERROR en predicción espacial:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (!is.null(occ_pred)) {
    # Insertar predicciones en el raster de referencia
    # NOTA: fit es la probabilidad de OCUPACIÓN (0 = no ocupado, 1 = totalmente ocupado)
    # El mapa muestra correctamente: valores ALTOS = mayor probabilidad de presencia
    r_pred <- rast(cov_stack[[1]])
    values(r_pred) <- NA_real_
    
    # Identificar celdas sin NA en todas las covariables
    celdas_validas <- complete.cases(pred_surface[, cov_names])
    values(r_pred)[celdas_validas] <- occ_pred$fit[celdas_validas]
    names(r_pred) <- "Probabilidad_Ocupacion"
    
    # Error estándar
    r_se <- rast(cov_stack[[1]])
    values(r_se) <- NA_real_
    values(r_se)[celdas_validas] <- occ_pred$se.fit[celdas_validas]
    names(r_se) <- "Error_Estandar"
    
    # Exportar rasters
    writeRaster(r_pred,
                file.path(path_spp, paste0("Mapa_Ocupacion_", slug, ".tif")),
                overwrite = TRUE)
    writeRaster(r_se,
                file.path(path_spp, paste0("Mapa_SE_Ocupacion_", slug, ".tif")),
                overwrite = TRUE)
    
    # Mapa clasificado (1–10)
    reclass_m <- matrix(
      c(0.0, 0.1, 1,  0.1, 0.2, 2,  0.2, 0.3, 3,
        0.3, 0.4, 4,  0.4, 0.5, 5,  0.5, 0.6, 6,
        0.6, 0.7, 7,  0.7, 0.8, 8,  0.8, 0.9, 9,
        0.9, 1.0, 10),
      ncol = 3, byrow = TRUE
    )
    r_clasificado <- classify(r_pred, reclass_m)
    writeRaster(r_clasificado,
                file.path(path_spp, paste0("Mapa_Ocupacion_Clasificado_", slug, ".tif")),
                overwrite = TRUE)
    
    # Gráfica del mapa
    # scale_fill_viridis_c: valores BAJOS = morado (ausente), valores ALTOS = amarillo (ocupado)
    mapa_plot <- ggplot() +
      geom_spatraster(data = r_pred) +
      scale_fill_viridis_c(
        option   = "viridis",
        na.value = "transparent",
        name     = "ψ (Prob. Ocupación)",
        limits   = c(0, 1),
        # direction = 1 asegura que 0=oscuro (ausente), 1=claro (ocupado)
        direction = 1
      ) +
      geom_sf(data = puntos_sf, color = "red", size = 0.8, alpha = 0.6) +
      labs(
        title    = paste("Mapa de Ocupación —", especie_objetivo),
        subtitle = "Modelo promediado (top 95% AICc weight) | Rojo = cámaras trampa",
        caption  = "0 = sin ocupación, 1 = máxima probabilidad de ocupación"
      ) +
      theme_minimal(base_size = 11)
    
    ggsave(
      filename = file.path(path_spp, paste0("Mapa_Ocupacion_", slug, ".png")),
      plot     = mapa_plot,
      width = 10, height = 8, dpi = 300
    )
    
    # Mapa de incertidumbre (SE)
    mapa_se_plot <- ggplot() +
      geom_spatraster(data = r_se) +
      scale_fill_viridis_c(
        option   = "magma",
        na.value = "transparent",
        name     = "Error estándar",
        direction = 1
      ) +
      labs(
        title    = paste("Incertidumbre (SE) —", especie_objetivo),
        subtitle = "Valores altos = mayor incertidumbre en la estimación"
      ) +
      theme_minimal(base_size = 11)
    
    ggsave(
      filename = file.path(path_spp, paste0("Mapa_SE_Ocupacion_", slug, ".png")),
      plot     = mapa_se_plot,
      width = 10, height = 8, dpi = 300
    )
    
  } else {
    r_pred <- NULL
    cat("  Mapa no generado por error en predicción.\n")
  }
  
  # --------------------------------------------------------
  # 3.11 INDICADORES FINALES
  # --------------------------------------------------------
  psi_sin_cov <- backTransform(fm0, type = "state")
  p_sin_cov   <- backTransform(fm0, type = "det")
  
  area_ocupada_pct <- if (!is.null(r_pred)) {
    global(r_pred, fun = "mean", na.rm = TRUE)[1, 1] * 100
  } else { NA_real_ }
  
  indicadores <- data.frame(
    Especie   = especie_objetivo,
    Indicador = c(
      "Probabilidad de Ocupación (ψ) — modelo nulo",
      "Probabilidad de Detección (p) — modelo nulo",
      "Porcentaje del área de estudio con ψ > 0.5",
      "N° modelos en promediado (95% weight)",
      "N° detecciones totales"
    ),
    Valor = c(
      round(psi_sin_cov@estimate, 3),
      round(p_sin_cov@estimate,   3),
      if (!is.null(r_pred)) round(global(r_pred > 0.5, fun = "sum", na.rm = TRUE)[1, 1] /
                                    global(!is.na(r_pred), fun = "sum", na.rm = TRUE)[1, 1] * 100, 2) else NA,
      length(occ_dredge_95),
      sum(y_final, na.rm = TRUE)
    ),
    Error_Estandar = c(
      round(SE(psi_sin_cov), 3),
      round(SE(p_sin_cov),   3),
      NA, NA, NA
    )
  )
  
  write.csv(indicadores,
            file.path(path_spp, paste0("Indicadores_Resumen_", slug, ".csv")),
            row.names = FALSE)
  
  cat("  Indicadores exportados:\n")
  print(indicadores)
  cat("  Resultados guardados en:", path_spp, "\n")
  
  return(invisible(list(
    modelo_avg   = occ_avg,
    dredge_table = tabla_dredge,
    indicadores  = indicadores,
    mapa         = r_pred
  )))
}

# ==========================================================
# 4. EJECUTAR ANÁLISIS PARA TODAS LAS ESPECIES
# ==========================================================
resultados <- list()

for (sp in especies_objetivo) {
  resultados[[sp]] <- tryCatch(
    analizar_especie(sp),
    error = function(e) {
      cat("\n  ERROR INESPERADO en", sp, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
}

cat("\n========================================================\n")
cat("  ANÁLISIS COMPLETADO PARA TODAS LAS ESPECIES\n")
cat("========================================================\n")

# ==========================================================
# 5. TABLA RESUMEN COMPARATIVA (TODAS LAS ESPECIES)
# ==========================================================
resumen_global <- lapply(nombres_con_resultado, function(sp) {
  res <- resultados[[sp]]
  if (is.null(res)) return(NULL)
  res$indicadores
}) |>
  bind_rows()

nombres_con_resultado <- names(Filter(Negate(is.null), resultados))

resumen_global <- lapply(nombres_con_resultado, function(sp) {
  resultados[[sp]]$indicadores
}) |> bind_rows()

write.csv(resumen_global,
          file.path(path_out, "Resumen_Global_Todas_Especies.csv"),
          row.names = FALSE)

cat("Tabla comparativa guardada en:", file.path(path_out, "Resumen_Global_Todas_Especies.csv"), "\n")