# ==============================================================================
# 00_correr_todo.R
# Corre los 4 scripts en orden. Es el unico archivo que necesitas ejecutar.
# Antes de correrlo, revisa que tu directorio de trabajo sea la raiz del
# repositorio (donde estan las carpetas data/ y scripts/).
# ==============================================================================

source("scripts/01_cargar_datos.R")
source("scripts/02_historia_deteccion.R")
source("scripts/03_modelo_ocupacion.R")
source("scripts/04_exportar_resultados.R")
