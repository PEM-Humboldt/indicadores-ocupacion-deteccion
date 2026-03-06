library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(tidyterra)

# ==========================================================
# 1. CONFIGURACIÓN DE DATOS Y ESPECIE
# ==========================================================

especie_objetivo <- "Nothocrax urumutum"
# Cambiamos las referencias a los objetos de fototrampeo
obs_data <- bioacustic_obs
dev_data <- bioacustic_dev
med_data <- bioacustic_media

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
         ocasion = as.numeric(floor(difftime(fecha, fecha_min, units = "days") / 7) + 1)) %>%
  group_by(deploymentID, ocasion) %>%
  summarise(presencia = 1, .groups = 'drop') %>%
  pivot_wider(names_from = ocasion, values_from = presencia, values_fill = 0)

# El error de "0, 82" se resuelve aquí:
y_final <- data.frame(deploymentID = site_covs_scaled$deploymentID) %>%
  left_join(y_matrix, by = "deploymentID") %>%
  select(-deploymentID) %>%
  as.matrix()

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
          "~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Nothocrax_urumutum.csv", 
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
         subtitle = "Modelo basado en datos de bioacustica",
         x = label_x, 
         y = "Probabilidad de Ocupación (psi)") +
    ylim(0, 1) # Obligamos a que el eje Y sea de 0 a 1 para ver la magnitud real
}

# --- EJECUTAR GRÁFICAS ---

# Guardar curva de Bosque
p_bosque <- plot_occu_tapir(fm1, "bosque", "Cobertura de Bosque (%)", "darkgreen")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_Bosque_Nothocrax_urumutum.png", 
       plot = p_bosque, width = 7, height = 5)

# Guardar curva de Ríos
p_rios <- plot_occu_tapir(fm2, "dist_rios", "Distancia a Ríos (m)", "royalblue")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_Rios_Nothocrax_urumutum.png", 
       plot = p_rios, width = 7, height = 5)

# Guardar curva de Vías
p_vias <- plot_occu_tapir(fm3, "dist_vias", "Distancia a Vías (m)", "firebrick")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_Vias_Nothocrax_urumutum.png", 
       plot = p_vias, width = 7, height = 5)


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
  scale_fill_viridis_c(option = "viridis", name = "Psi (Nothocrax_urumutum)", limits = c(0, 1)) +
  geom_sf(data = puntos_sf, color = "red", size = 0.5) +
  labs(title = paste("Ocupación Proyectada:", especie_objetivo)) +
  theme_minimal()


# 1. Creamos el objeto del mapa
mapa_final <- ggplot() +
  geom_spatraster(data = mapa_ocupacion) +
  scale_fill_viridis_c(option = "magma", na.value = "transparent",
                       name = "Psi (Nothocrax_urumutum)", limits = c(0, 1)) +
  geom_sf(data = puntos_sf, color = "white", size = 1, alpha = 0.5) +
  # ... (aquí van tus labels, norte, escala y theme que ajustamos antes) ...
  labs(title = paste("Mapa de Ocupación:", especie_objetivo)) +
  theme_minimal()

# 2. Exportamos especificando el objeto 'plot = mapa_final'
ggsave(filename = "~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Mapa_Final_Nothocrax_urumutum.png", 
       plot = mapa_final, 
       width = 10, 
       height = 8, 
       dpi = 300)


# 1. Guardar el Raster de Probabilidad (TIFF)
# Este es el archivo clave para que lo abras en QGIS
writeRaster(mapa_ocupacion, 
            "~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Mapa_Ocupacion_Nothocrax_urumutum.tif", 
            overwrite = TRUE)
