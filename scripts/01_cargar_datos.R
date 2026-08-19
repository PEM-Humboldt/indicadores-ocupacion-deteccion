# ==============================================================================
# 01_cargar_datos.R
# Carga las observaciones, las covariables de sitio y el grid de prediccion.
#
# UNICO PUNTO QUE DEBES EDITAR: la ruta del proyecto (path_base).
# ==============================================================================

library(dplyr)
library(readr)

# --- 1. Rutas -----------------------------------------------------------
path_base <- "."

path_obs  <- file.path(path_base, "data/observaciones.csv")
path_cov  <- file.path(path_base, "data/sitios_covariables.csv")
path_grid <- file.path(path_base, "data/grid_prediccion.csv")
path_out  <- file.path(path_base, "resultados")

if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

# --- 2. Cargar tablas -----------------------------------------------------
observaciones <- read_csv(path_obs, show_col_types = FALSE)
covariables_sitio <- read_csv(path_cov, show_col_types = FALSE)
grid_prediccion <- read_csv(path_grid, show_col_types = FALSE)

cat("=== DATOS CARGADOS ===\n")
cat("Sitios de muestreo :", n_distinct(observaciones$sitio), "\n")
cat("Visitas por sitio  :", max(observaciones$visita), "\n")
cat("Especies           :", paste(unique(observaciones$especie), collapse = ", "), "\n")
cat("Covariables        : pct_bosque, dist_rio_m, dist_via_m\n")
cat("Celdas en el grid de prediccion:", nrow(grid_prediccion), "\n")
