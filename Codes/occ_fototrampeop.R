library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(tidyterra)
library(data.table)

# ==========================================================
# 1. CONFIGURACIÓN DE DATOS Y ESPECIE
# ==========================================================

# Cambiamos las referencias a los objetos de fototrampeo
path_data <- "~/Desktop/FPVA/Data/Fototrampeo/I2D_FPVA_Fototrampeo_20260219.xlsx"
path_out  <- "~/Desktop/FPVA/Resultados/Fototrampeo/"
if(!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

obs_data  <- read.xlsx(path_data, sheet = "Observations")
med_data <- read.xlsx(path_data, sheet = "Media")
dev_data  <- read.xlsx(path_data, sheet = "Deployment")

especie_objetivo <- "Tapirus terrestris"
# ==========================================================
# 2. PREPARACIÓN DE COVARIABLES (Limpieza y Estandarización)
# ==========================================================
path_base <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/Co-variables Ocupación/"

cov_stack <- c(
  rast(paste0(path_base, "BosqueNoBosque__90mts/agreg_bosque90_2.tif")),
  rast(paste0(path_base, "EucDistance_Drenajes90mts/RiosEudis90.tif")),
  rast(paste0(path_base, "EucliDistance_vias90mts/ViasEudis90.tif"))
)
names(cov_stack) <- c("bosque", "dist_rios", "dist_vias")

# Extraer usando los puntos de FOTOTRAMPEO
puntos_sf <- dev_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(cov_stack))

site_covs_raw <- terra::extract(cov_stack, puntos_sf)
site_covs_raw$deploymentID <- dev_data$deploymentID

# Guardar parámetros
params <- site_covs_raw %>%
  summarise(across(c(bosque, dist_rios, dist_vias), 
                   list(mu = ~mean(.x, na.rm=T), sd = ~sd(.x, na.rm=T))))

site_covs_scaled <- site_covs_raw %>%
  mutate(across(c(bosque, dist_rios, dist_vias), scale)) %>%
  filter(complete.cases(bosque, dist_rios, dist_vias))

# ==========================================================
# 3. HISTORIA DE DETECCIÓN (Tapirus terrestris)
# ==========================================================
# Unimos observaciones con media para tener la estampa de tiempo
# IMPORTANTE: Asegúrate que mediaID sea la llave de unión
obs_con_tiempo <- obs_data %>%
  left_join(med_data %>% select(mediaID, timestamp), by = "mediaID") %>%
  mutate(time_obj = ymd_hms(timestamp)) %>%
  filter(!is.na(time_obj))

fecha_min <- min(as.Date(obs_con_tiempo$time_obj), na.rm = TRUE)

y_matrix <- obs_con_tiempo %>%
  filter(scientificName == especie_objetivo) %>%
  mutate(fecha = as.Date(time_obj),
         ocasion = as.numeric(floor(difftime(fecha, fecha_min, units = "days") / 15) + 1)) %>%
  group_by(deploymentID, ocasion) %>%
  summarise(presencia = 1, .groups = 'drop') %>%
  pivot_wider(names_from = ocasion, values_from = presencia, values_fill = 0)

# El error de "0, 82" se resuelve aquí:
y_final <- data.frame(deploymentID = site_covs_scaled$deploymentID) %>%
  left_join(y_matrix, by = "deploymentID") %>%
  select(-deploymentID) %>%
  as.matrix()
# revisar esta linea
y_final[is.na(y_final)] <- 0 

# ==========================================================
# 4. MODELADO unmarked
# ==========================================================
site_covs_unmarked <- site_covs_scaled %>% select(bosque, dist_rios, dist_vias)
umf <- unmarkedFrameOccu(y = y_final, siteCovs = site_covs_unmarked)

# Ajuste de modelos
fm0 <- occu(~1 ~1, data = umf)
fm1 <- occu(~1 ~bosque, data = umf)
fm2 <- occu(~1 ~dist_rios, data = umf)
fm3 <- occu(~1 ~dist_vias, data = umf)
fm_global <- occu(~1 ~bosque + dist_rios + dist_vias, data = umf)

model_list <- fitList('Null'=fm0, 'Bosque'=fm1, 'Rios'=fm2, 'Vias'=fm3, 'Global'=fm_global)
print(modSel(model_list))

# Extraer la tabla de selección como data.frame correctamente
tabla_aic <- as(modSel(model_list), "data.frame")

# Exportar ahora sí sin errores
write.csv(tabla_aic, 
          "~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Tabla_AIC_Modelos_Dasypus_novemcinctus.csv", 
          row.names = FALSE)
# ==========================================================
# 4. GRÁFICAS DE RESPUESTA PARA TAPIRUS TERRESTRIS
# ==========================================================

# Función maestra ajustada para fototrampeo
plot_occu_tapir <- function(model, var_name, label_x, color_fill) {
  
  # 1. Recuperar parámetros originales de la tabla de fototrampeo
  mu <- params[[paste0(var_name, "_mu")]]
  sd <- params[[paste0(var_name, "_sd")]]
  
  # 2. Rango real (de min a max observado en los sitios de cámaras)
  min_val <- min(site_covs_raw[[var_name]], na.rm=T)
  max_val <- max(site_covs_raw[[var_name]], na.rm=T)
  raw_seq <- seq(min_val, max_val, length = 100)
  
  # 3. Escalar para que el modelo lo entienda
  scaled_seq <- (raw_seq - mu) / sd
  
  # Crear dataframe de predicción
  nd <- data.frame(scaled_seq)
  names(nd) <- var_name
  
  # Predecir ocupación (psi)
  pred <- predict(model, type = "state", newdata = nd, appendData = TRUE)
  
  # Ajustar eje X si es bosque para mostrar porcentaje (0-100)
  x_display <- if(var_name == "bosque") raw_seq * 100 else raw_seq
  
  # 4. Generar el gráfico
  ggplot(pred, aes(x = x_display, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = color_fill) +
    geom_line(color = color_fill, size = 1.2) +
    theme_minimal() +
    labs(title = paste("Respuesta de", especie_objetivo, "a", label_x),
         subtitle = "Modelo basado en datos de fototrampeo",
         x = label_x, 
         y = "Probabilidad de Ocupación (psi)") +
    ylim(0, 1) # Obligamos a que el eje Y sea de 0 a 1 para ver la magnitud real
}

# --- EJECUTAR GRÁFICAS ---

# Guardar curva de Bosque
p_bosque <- plot_occu_tapir(fm1, "bosque", "Cobertura de Bosque (%)", "darkgreen")
ggsave("~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Curva_Bosque_Dasypus_novemcinctus.png", 
       plot = p_bosque, width = 7, height = 5)

# Guardar curva de Ríos
p_rios <- plot_occu_tapir(fm2, "dist_rios", "Distancia a Ríos (m)", "royalblue")
ggsave("~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Curva_Rios_Dasypus_novemcinctus.png", 
       plot = p_rios, width = 7, height = 5)

# Guardar curva de Vías
p_vias <- plot_occu_tapir(fm3, "dist_vias", "Distancia a Vías (m)", "firebrick")
ggsave("~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Curva_Vias_Dasypus_novemcinctus.png", 
       plot = p_vias, width = 7, height = 5)

# ==========================================================
# 4.1 PREPARACIÓN DE COVARIABLES DE DETECCIÓN (obsCovs)
# ==========================================================

# 4.1.1 Matriz de Esfuerzo (Días activos por ocasión de 7 días)
# Esto corrige el sesgo de cámaras que se apagaron antes
esfuerzo_matrix <- obs_con_tiempo %>%
  group_by(deploymentID) %>%
  summarise(inicio_site = min(as.Date(time_obj)),
            fin_site = max(as.Date(time_obj))) %>%
  right_join(data.frame(deploymentID = site_covs_scaled$deploymentID), by = "deploymentID")

# Creamos una matriz de esfuerzo (días por semana)
# (Para este ejemplo simplificado, asumiremos 7 días si la cámara estaba activa, 
# pero lo ideal es cruzarlo con tu tabla 'dev' para ver fechas exactas de falla)
n_ocasiones <- ncol(y_final)
effort_mat <- matrix(7, nrow = nrow(y_final), ncol = n_ocasiones)

# 4.2. Hora de detección (Promedio de hora por sitio)
# La usamos como siteCov porque describe el comportamiento en ese lugar
hora_sitio <- obs_con_tiempo %>%
  filter(scientificName == especie_objetivo) %>%
  mutate(hora = hour(time_obj)) %>%
  group_by(deploymentID) %>%
  summarise(hora_promedio = mean(hora, na.rm = TRUE)) %>%
  right_join(data.frame(deploymentID = site_covs_scaled$deploymentID), by = "deploymentID") %>%
  mutate(hora_promedio = replace_na(hora_promedio, mean(hora_promedio, na.rm = TRUE))) # Imputar media global a sitios sin detección

site_covs_unmarked <- site_covs_scaled %>% 
  mutate(hora_actividad = scale(hora_sitio$hora_promedio)) %>%
  select(bosque, dist_rios, dist_vias, hora_actividad)

# ==========================================================
# 4.2 MODELADO CON DETECCIÓN VARIABLE
# ==========================================================

# Creamos el unmarkedFrame con siteCovs y obsCovs
umf <- unmarkedFrameOccu(
  y = y_final, 
  siteCovs = site_covs_unmarked,
  obsCovs = list(esfuerzo = effort_mat)
)

# Modelos comparando efectos en Ocupación (psi) y Detección (p)
# Formato: ~ deteccion ~ ocupacion

# 1. Modelo Nulo (p y psi constantes)
fm0 <- occu(~1 ~1, data = umf)

# 2. ¿Afecta el bosque la ocupación? (p constante)
fm1 <- occu(~1 ~bosque, data = umf)

# 3. ¿Afecta el bosque la detección? (psi constante)
fm2 <- occu(~bosque ~1, data = umf)

# 4. ¿Afecta la hora de actividad la detección?
fm3 <- occu(~hora_actividad ~1, data = umf)

# 5. Modelo Global (Efectos en ambos parámetros)
fm_global <- occu(~bosque + hora_actividad ~bosque + dist_rios + dist_vias, data = umf)

# Selección de modelos
model_list <- fitList(
  'Nulo' = fm0, 
  'Psi(Bosque)' = fm1, 
  'P(Bosque)' = fm2, 
  'P(Hora)' = fm3, 
  'Global' = fm_global
)

tab_aic <- as(modSel(model_list), "data.frame")
print(tab_aic)

# ==========================================================
# 4.3. GRÁFICA DE RESPUESTA DE DETECCIÓN (p)
# ==========================================================

plot_det_prob <- function(model, var_name, label_x) {
  # Crear secuencia de predicción para la detección (type = "det")
  raw_seq <- seq(min(site_covs_unmarked[[var_name]]), max(site_covs_unmarked[[var_name]]), length = 100)
  nd <- data.frame(raw_seq)
  names(nd) <- var_name
  
  pred <- predict(model, type = "det", newdata = nd, appendData = TRUE)
  
  ggplot(pred, aes(x = raw_seq, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "purple") +
    geom_line(color = "purple", size = 1.2) +
    theme_minimal() +
    labs(title = paste("Probabilidad de Detección (p) de", especie_objetivo),
         x = label_x, y = "Probabilidad de Detección (p)") +
    ylim(0, 1)
}

# Ejemplo: Graficar cómo afecta el bosque a la detección
if ("P(Bosque)" %in% names(model_list)) {
  p_det_bosque <- plot_det_prob(fm2, "bosque", "Cobertura de Bosque (Escalado)")
  ggsave("~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Deteccion_Bosque.png", p_det_bosque)
}

####################
# Graficas detección
#######################

# Asegúrate de que los nombres en fitList coincidan exactamente aquí
# Si no quieres complicaciones con el 'if', simplemente corre la línea directamente:

# 1. Definir la función mejorada (des-escalando el eje X)
plot_det_prob_real <- function(model, var_name, label_x) {
  
  # Recuperar parámetros para des-escalar (mu y sd)
  mu <- params[[paste0(var_name, "_mu")]]
  sd <- params[[paste0(var_name, "_sd")]]
  
  # Secuencia en valores escalados (lo que entiende el modelo)
  scaled_seq <- seq(min(site_covs_unmarked[[var_name]]), 
                    max(site_covs_unmarked[[var_name]]), length = 100)
  
  nd <- data.frame(scaled_seq)
  names(nd) <- var_name
  
  # Predicción de detección (type = "det")
  pred <- predict(model, type = "det", newdata = nd, appendData = TRUE)
  
  # Volver a valores reales para el gráfico
  x_real <- (scaled_seq * sd) + mu
  if(var_name == "bosque") x_real <- x_real * 100 # Convertir a %
  
  ggplot(pred, aes(x = x_real, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "purple") +
    geom_line(color = "purple", size = 1.2) +
    theme_minimal() +
    labs(title = paste("Detección (p) de", especie_objetivo),
         subtitle = paste("Efecto de:", label_x),
         x = label_x, y = "Probabilidad de Detección (p)") +
    ylim(0, 1)
}

# 2. Ejecutar directamente (sin el IF para evitar el error de objeto no encontrado)
p_det_bosque <- plot_det_prob_real(fm2, "bosque", "Cobertura de Bosque (%)")

# 3. Guardar
ggsave("~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Deteccion_Bosque.png", 
       plot = p_det_bosque, width = 7, height = 5)

# Ocupación (tiene 4 coeficientes: Int, bosque, rios, vias)
psi_media <- backTransform(linearComb(fm_global, coefficients = c(1,0,0,0), type = "state"))@estimate

# Detección (tiene 3 coeficientes: Int, bosque, hora_actividad)
# Ajustamos a c(1,0,0) para que coincida con el número de variables en 'p'
p_media <- backTransform(linearComb(fm_global, coefficients = c(1,0,0), type = "det"))@estimate

cat("Ocupación media (psi):", psi_media, "\n")
cat("Detección media (p):", p_media, "\n")

####
#
# 1. Asegurémonos de que 'params' tenga la información de 'hora_actividad'
params$hora_actividad_mu <- mean(hora_sitio$hora_promedio, na.rm = TRUE)
params$hora_actividad_sd <- sd(hora_sitio$hora_promedio, na.rm = TRUE)

# 2. Función corregida
plot_det_v3 <- function(model, var_name, label_x, color_line) {
  
  # Rango de la variable escalada (basado en lo que el modelo vio)
  vals_escalados <- seq(min(site_covs_unmarked[[var_name]]), 
                        max(site_covs_unmarked[[var_name]]), length = 100)
  
  # Crear el newdata con TODAS las variables que tiene la parte 'p' del modelo
  # En fm_global, la detección es: ~bosque + hora_actividad
  if(var_name == "bosque") {
    nd <- data.frame(bosque = vals_escalados, hora_actividad = 0)
  } else {
    nd <- data.frame(bosque = 0, hora_actividad = vals_escalados)
  }
  
  # Predecir
  pred <- predict(model, type = "det", newdata = nd, appendData = TRUE)
  
  # Recuperar mu y sd para des-escalar
  mu <- params[[paste0(var_name, "_mu")]]
  sd <- params[[paste0(var_name, "_sd")]]
  
  # Añadir el eje X real al dataframe de resultados para que ggplot no se pierda
  pred$x_real <- (vals_escalados * sd) + mu
  
  # Ajuste estético para bosque (%)
  if(var_name == "bosque") pred$x_real <- pred$x_real * 100
  
  # Graficar
  ggplot(pred, aes(x = x_real, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = color_line) +
    geom_line(color = color_line, size = 1.2) +
    theme_minimal() +
    labs(title = paste("Detección (p) de", especie_objetivo),
         subtitle = paste("Efecto de:", label_x),
         x = label_x, y = "Probabilidad de Detección") +
    ylim(0, 1)
}

# --- EJECUTAR ---
p_det_hora <- plot_det_v3(fm_global, "hora_actividad", "Hora Promedio de Actividad (Decimal)", "darkorange")
p_det_bosque <- plot_det_v3(fm_global, "bosque", "Cobertura de Bosque (%)", "darkgreen")

# Visualizar con patchwork
library(patchwork)
p_det_hora / p_det_bosque

# ==========================================================
# 5. MAPA DE OCUPACIÓN (Proyección espacial)
# ==========================================================
# Seleccionamos el mejor modelo (supongamos el global)
modelo_ganador <- fm_global 
beta <- coef(modelo_ganador, type = "state")

# Escalar el stack para predicción
env_stack_scaled <- cov_stack
env_stack_scaled$bosque    <- (cov_stack$bosque - params$bosque_mu) / params$bosque_sd
env_stack_scaled$dist_rios <- (cov_stack$dist_rios - params$dist_rios_mu) / params$dist_rios_sd
env_stack_scaled$dist_vias <- (cov_stack$dist_vias - params$dist_vias_mu) / params$dist_vias_sd

# Cálculo de probabilidad
logit_psi <- beta[1] + 
  beta[2] * env_stack_scaled$bosque + 
  beta[3] * env_stack_scaled$dist_rios + 
  beta[4] * env_stack_scaled$dist_vias

mapa_ocupacion <- exp(logit_psi) / (1 + exp(logit_psi))
names(mapa_ocupacion) <- "Probabilidad_Ocupacion"

# Graficar
ggplot() +
  geom_spatraster(data = mapa_ocupacion) +
  scale_fill_viridis_c(option = "viridis", name = "Psi (Dasypus_novemcinctus)", limits = c(0, 1)) +
  geom_sf(data = puntos_sf, color = "red", size = 0.5) +
  labs(title = paste("Ocupación Proyectada:", especie_objetivo)) +
  theme_minimal()


# 1. Creamos el objeto del mapa
mapa_final <- ggplot() +
  geom_spatraster(data = mapa_ocupacion) +
  scale_fill_viridis_c(option = "magma", na.value = "transparent",
                       name = "Psi (Dasypus_novemcinctus)", limits = c(0, 1)) +
  geom_sf(data = puntos_sf, color = "white", size = 1, alpha = 0.5) +
  # ... (aquí van tus labels, norte, escala y theme que ajustamos antes) ...
  labs(title = paste("Mapa de Ocupación:", especie_objetivo)) +
  theme_minimal()

# 2. Exportamos especificando el objeto 'plot = mapa_final'
ggsave(filename = "~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Mapa_Final_Dasypus_novemcinctus.png", 
       plot = mapa_final, 
       width = 10, 
       height = 8, 
       dpi = 300)


# 1. Guardar el Raster de Probabilidad (TIFF)
# Este es el archivo clave para que lo abras en QGIS
writeRaster(mapa_ocupacion, 
            "~/Desktop/FPVA/Resultados/Fototrampeo/Dasypus_novemcinctus/Mapa_Ocupacion_Dasypus_novemcinctus.tif", 
            overwrite = TRUE)
