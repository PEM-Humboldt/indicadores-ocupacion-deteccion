# ==============================================================================
# 03_modelo_ocupacion.R
# Indicador 1: probabilidad de ocupacion (psi) sin covariables
# Indicador 2: probabilidad de deteccion (p) sin covariables
# Indicador 3: porcentaje del area ocupada (proyectando el modelo con
#              covariables sobre el grid de prediccion)
#
# Covariables (siempre las mismas 3, estandarizadas a media 0 / sd 1):
#   pct_bosque, dist_rio_m, dist_via_m
#
# Requiere haber corrido antes: 01_cargar_datos.R y 02_historia_deteccion.R
# ==============================================================================

library(dplyr)
library(unmarked)

UMBRAL_PSI <- 0.5  # umbral para el % de area ocupada (indicador 3)

# --- 1. Estandarizar covariables (media 0, sd 1) ---------------------------
estandarizar <- function(x, media, sd) (x - media) / sd

medias <- sapply(covariables_sitio[, c("pct_bosque", "dist_rio_m", "dist_via_m")], mean)
sds    <- sapply(covariables_sitio[, c("pct_bosque", "dist_rio_m", "dist_via_m")], sd)

cov_std <- covariables_sitio |>
  transmute(
    pct_bosque = estandarizar(pct_bosque, medias["pct_bosque"], sds["pct_bosque"]),
    dist_rio_m = estandarizar(dist_rio_m, medias["dist_rio_m"], sds["dist_rio_m"]),
    dist_via_m = estandarizar(dist_via_m, medias["dist_via_m"], sds["dist_via_m"])
  )

# --- 2. Armar el objeto unmarkedFrameOccu -----------------------------------
y_mat <- as.matrix(historia_deteccion[, -1])

umf <- unmarkedFrameOccu(y = y_mat, siteCovs = cov_std)

# --- 3. Modelo nulo (sin covariables) — indicadores 1 y 2 -------------------
modelo_nulo <- occu(~1 ~1, data = umf)

psi_nulo <- round(plogis(coef(modelo_nulo)["psi(Int)"]), 3)
p_nulo   <- round(plogis(coef(modelo_nulo)["p(Int)"]), 3)

cat("=== MODELO NULO (sin covariables) ===\n")
cat("psi (ocupacion) :", psi_nulo, "\n")
cat("p (deteccion)   :", p_nulo, "\n")

# --- 4. Modelo con las 3 covariables en psi, deteccion constante -----------
modelo_cov <- occu(~1 ~ pct_bosque + dist_rio_m + dist_via_m, data = umf)

# --- 5. Seleccion de modelo por AIC -----------------------------------------
tab_aic <- data.frame(
  modelo = c("Nulo (psi~1, p~1)", "Covariables (psi~bosque+rio+via, p~1)"),
  AIC    = c(AIC(modelo_nulo), AIC(modelo_cov))
) |> arrange(AIC)

cat("\n=== SELECCION DE MODELO (AIC) ===\n")
print(tab_aic)

modelo_final <- if (tab_aic$modelo[1] == "Nulo (psi~1, p~1)") modelo_nulo else modelo_cov
cat("\nModelo seleccionado:", tab_aic$modelo[1], "\n")

# --- 6. Indicador 3: % de area ocupada, proyectando sobre el grid ----------
grid_std <- grid_prediccion |>
  transmute(
    id_celda   = id_celda,
    pct_bosque = estandarizar(pct_bosque, medias["pct_bosque"], sds["pct_bosque"]),
    dist_rio_m = estandarizar(dist_rio_m, medias["dist_rio_m"], sds["dist_rio_m"]),
    dist_via_m = estandarizar(dist_via_m, medias["dist_via_m"], sds["dist_via_m"])
  )

pred_grid <- predict(modelo_cov, type = "state", newdata = grid_std)

grid_prediccion$psi_predicho <- round(pred_grid$Predicted, 3)
grid_prediccion$psi_SE       <- round(pred_grid$SE, 3)

pct_area_ocupada <- round(mean(grid_prediccion$psi_predicho >= UMBRAL_PSI) * 100, 1)

cat("\n=== INDICADOR 3: % DEL AREA OCUPADA ===\n")
cat("Umbral psi >=", UMBRAL_PSI, "\n")
cat("Celdas del grid con psi >= umbral:",
    sum(grid_prediccion$psi_predicho >= UMBRAL_PSI), "de", nrow(grid_prediccion), "\n")
cat("Porcentaje de area ocupada:", pct_area_ocupada, "%\n")

# --- 7. Tabla resumen de los 3 indicadores ----------------------------------
tab_indicadores <- data.frame(
  especie = especie_objetivo,
  psi_nulo = psi_nulo,
  p_nulo = p_nulo,
  modelo_seleccionado = tab_aic$modelo[1],
  pct_area_ocupada = pct_area_ocupada,
  umbral_psi = UMBRAL_PSI
)

cat("\n=== RESUMEN ===\n")
print(tab_indicadores)
