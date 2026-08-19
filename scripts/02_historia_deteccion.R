# ==============================================================================
# 02_historia_deteccion.R
# Construye la historia de detección (matriz sitio x visita, con 0/1) que
# necesita el paquete `unmarked` para ajustar el modelo de ocupación.
#
# Requiere haber corrido antes: 01_cargar_datos.R
# ==============================================================================

library(dplyr)
library(tidyr)

# Edita este nombre si quieres modelar otra especie de tu tabla de datos
especie_objetivo <- "Dicotyles tajacu"

historia_deteccion <- observaciones |>
  filter(especie == especie_objetivo) |>
  select(sitio, visita, detectado) |>
  pivot_wider(names_from = visita, values_from = detectado,
              names_prefix = "visita_") |>
  arrange(sitio)

# Aseguramos que las covariables queden en el mismo orden de sitios que
# la historia de deteccion (paso clave para que unmarked no las cruce mal)
covariables_sitio <- covariables_sitio |>
  filter(sitio %in% historia_deteccion$sitio) |>
  arrange(sitio)

stopifnot(identical(historia_deteccion$sitio, covariables_sitio$sitio))

cat("=== HISTORIA DE DETECCION —", especie_objetivo, "===\n")
print(historia_deteccion)

cat("\nSitios con al menos una deteccion:",
    sum(rowSums(historia_deteccion[,-1], na.rm = TRUE) > 0),
    "de", nrow(historia_deteccion), "\n")
