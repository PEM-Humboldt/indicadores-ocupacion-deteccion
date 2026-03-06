library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(tidyterra)
library(openxlsx)

# ==========================================================
# 1. CONFIGURACIÓN DE DATOS Y ESPECIE
# ==========================================================

path_data <- "~/Desktop/FPVA/Data/Puntos_conteo/I2D-BIO_2025_023.xlsx"
path_out  <- "~/Desktop/FPVA/Resultados/Aves//"
if(!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

obs_data  <- read.xlsx(path_data, sheet = "Registros")   # Tu primera tabla (registros)
med_data <- read.xlsx(path_data, sheet = "Evento")   # Tu segunda tabla (eventos/metadatos)

especie_objetivo <- "Cyanocorax violaceus" # Ejemplo basado en tus datos

# ==========================================================
# 2. PREPARACIÓN DE COVARIABLES DE SITIO (Rasters)
# ==========================================================
path_base <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/Co-variables Ocupación/"

cov_stack <- c(
  rast(paste0(path_base, "BosqueNoBosque__90mts/agreg_bosque90_2.tif")),
  rast(paste0(path_base, "EucDistance_Drenajes90mts/RiosEudis90.tif")),
  rast(paste0(path_base, "EucliDistance_vias90mts/ViasEudis90.tif"))
)
names(cov_stack) <- c("bosque", "dist_rios", "dist_vias")

# Extraer usando los puntos de AVES (parentEventID representa el sitio único)
sitios_coords <- med_data %>%
  group_by(parentEventID) %>%
  summarise(decimalLatitude = first(decimalLatitude),
            decimalLongitude = first(decimalLongitude))

puntos_sf <- sitios_coords %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs(cov_stack))

site_covs_raw <- terra::extract(cov_stack, puntos_sf)
site_covs_raw$parentEventID <- sitios_coords$parentEventID

# Guardar parámetros para des-escalar
params <- site_covs_raw %>%
  summarise(across(c(bosque, dist_rios, dist_vias), 
                   list(mu = ~mean(.x, na.rm=T), sd = ~sd(.x, na.rm=T))))

site_covs_scaled <- site_covs_raw %>%
  mutate(across(c(bosque, dist_rios, dist_vias), scale)) %>%
  filter(complete.cases(bosque, dist_rios, dist_vias))

# ==========================================================
# 3. HISTORIA DE DETECCIÓN (Basada en Visitas)
# ==========================================================
# En aves, cada eventID es una réplica (visita)
y_matrix <- obs_data %>%
  filter(scientificName == especie_objetivo) %>%
  # Creamos un número de visita por sitio
  group_by(parentEventID) %>%
  mutate(ocasion = as.numeric(as.factor(eventID))) %>%
  group_by(parentEventID, ocasion) %>%
  summarise(presencia = 1, .groups = 'drop') %>%
  pivot_wider(names_from = ocasion, values_from = presencia, values_fill = 0)

# Asegurar match con sitios que tienen covariables
y_final <- data.frame(parentEventID = site_covs_scaled$parentEventID) %>%
  left_join(y_matrix, by = "parentEventID") %>%
  select(-parentEventID) %>%
  as.matrix()

y_final[is.na(y_final)] <- 0 

# ==========================================================
# 4. COVARIABLES DE DETECCIÓN (Hora de inicio del punto)
# ==========================================================
# Extraemos la hora de inicio de eventTime (ej: "06:51/07:01" -> 6.85)
hora_deteccion <- med_data %>%
  mutate(hora_str = sub("/.*", "", eventTime),
         hora_obj = hm(hora_str),
         hora_dec = hour(hora_obj) + minute(hora_obj)/60) %>%
  group_by(parentEventID) %>%
  summarise(hora_mu = mean(hora_dec, na.rm=T))

params$hora_actividad_mu <- mean(hora_deteccion$hora_mu, na.rm=T)
params$hora_actividad_sd <- sd(hora_deteccion$hora_mu, na.rm=T)

site_covs_unmarked <- site_covs_scaled %>%
  left_join(hora_deteccion, by = "parentEventID") %>%
  mutate(hora_actividad = scale(hora_mu)) %>%
  select(bosque, dist_rios, dist_vias, hora_actividad)

# ==========================================================
# 5. MODELADO UNMARKED
# ==========================================================
umf <- unmarkedFrameOccu(y = y_final, siteCovs = site_covs_unmarked)

fm_global <- occu(~hora_actividad + bosque ~bosque + dist_rios + dist_vias, data = umf)

# Cálculos de probabilidad media
psi_media <- backTransform(linearComb(fm_global, coefficients = c(1,0,0,0), type = "state"))@estimate
p_media <- backTransform(linearComb(fm_global, coefficients = c(1,0,0), type = "det"))@estimate

cat("Ocupación media (psi):", psi_media, "\n")
cat("Detección media (p):", p_media, "\n")

# ==========================================================
# 6. MAPA Y EXPORTACIÓN
# ==========================================================
beta <- coef(fm_global, type = "state")
env_stack_scaled <- cov_stack
env_stack_scaled$bosque    <- (cov_stack$bosque - params$bosque_mu) / params$bosque_sd
env_stack_scaled$dist_rios <- (cov_stack$dist_rios - params$dist_rios_mu) / params$dist_rios_sd
env_stack_scaled$dist_vias <- (cov_stack$dist_vias - params$dist_vias_mu) / params$dist_vias_sd

logit_psi <- beta[1] + beta[2]*env_stack_scaled$bosque + 
  beta[3]*env_stack_scaled$dist_rios + beta[4]*env_stack_scaled$dist_vias
mapa_ocupacion <- exp(logit_psi) / (1 + exp(logit_psi))

# Guardar resultados
writeRaster(mapa_ocupacion, paste0("~/Desktop/FPVA/Resultados/Aves/", especie_objetivo, "_Ocupacion.tif"), overwrite=T)

# Graficar mapa final
ggplot() +
  geom_spatraster(data = mapa_ocupacion) +
  scale_fill_viridis_c(option = "magma", name = "Psi") +
  geom_sf(data = puntos_sf, color = "cyan", size = 1) +
  labs(title = paste("Mapa Ocupación Ave:", especie_objetivo)) +
  theme_minimal()

# ==========================================================
# 7. GRÁFICAS DE RESPUESTA: OCUPACIÓN (psi)
# ==========================================================

plot_occu_aves <- function(model, var_name, label_x, color_fill) {
  # Recuperar parámetros para des-escalar
  mu <- params[[paste0(var_name, "_mu")]]
  sd <- params[[paste0(var_name, "_sd")]]
  
  # Secuencia escalada
  scaled_seq <- seq(min(site_covs_unmarked[[var_name]]), 
                    max(site_covs_unmarked[[var_name]]), length = 100)
  
  # Preparar newdata (poniendo las otras variables en 0)
  if(var_name == "bosque") {
    nd <- data.frame(bosque = scaled_seq, dist_rios = 0, dist_vias = 0)
  } else if(var_name == "dist_rios") {
    nd <- data.frame(bosque = 0, dist_rios = scaled_seq, dist_vias = 0)
  } else {
    nd <- data.frame(bosque = 0, dist_rios = 0, dist_vias = scaled_seq)
  }
  
  # Predecir ocupación
  pred <- predict(model, type = "state", newdata = nd, appendData = TRUE)
  
  # Eje X real
  pred$x_real <- (scaled_seq * sd) + mu
  if(var_name == "bosque") pred$x_real <- pred$x_real * 100
  
  ggplot(pred, aes(x = x_real, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = color_fill) +
    geom_line(color = color_fill, size = 1.2) +
    theme_minimal() +
    labs(title = paste("Ocupación de", especie_objetivo),
         x = label_x, y = "Probabilidad de Ocupación (psi)") +
    ylim(0, 1)
}

# ==========================================================
# 8. GRÁFICAS DE RESPUESTA: DETECCIÓN (p)
# ==========================================================

plot_det_aves <- function(model, var_name, label_x, color_line) {
  scaled_seq <- seq(min(site_covs_unmarked[[var_name]]), 
                    max(site_covs_unmarked[[var_name]]), length = 100)
  
  # El modelo de detección es: ~hora_actividad + bosque
  if(var_name == "hora_actividad") {
    nd <- data.frame(hora_actividad = scaled_seq, bosque = 0)
  } else {
    nd <- data.frame(hora_actividad = 0, bosque = scaled_seq)
  }
  
  pred <- predict(model, type = "det", newdata = nd, appendData = TRUE)
  
  mu <- params[[paste0(var_name, "_mu")]]
  sd <- params[[paste0(var_name, "_sd")]]
  pred$x_real <- (scaled_seq * sd) + mu
  
  ggplot(pred, aes(x = x_real, y = Predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = color_line) +
    geom_line(color = color_line, size = 1.2) +
    theme_minimal() +
    labs(title = paste("Detección de", especie_objetivo),
         x = label_x, y = "Probabilidad de Detección (p)") +
    ylim(0, 1)
}

# ==========================================================
# 9. EJECUCIÓN Y EXPORTACIÓN DE PANELES
# ==========================================================
library(patchwork)

# Generar curvas de Ocupación
p1 <- plot_occu_aves(fm_global, "bosque", "Cobertura de Bosque (%)", "darkgreen")
p2 <- plot_occu_aves(fm_global, "dist_rios", "Distancia a Ríos (m)", "royalblue")
p3 <- plot_occu_aves(fm_global, "dist_vias", "Distancia a Vías (m)", "firebrick")

panel_psi <- p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave(paste0("~/Desktop/FPVA/Resultados/Aves/Respuesta_Ocupacion_", especie_objetivo, ".png"), 
       panel_psi, width = 15, height = 5)

# Generar curvas de Detección
p4 <- plot_det_aves(fm_global, "hora_actividad", "Hora de Inicio (Censo)", "darkorange")
p5 <- plot_det_aves(fm_global, "bosque", "Cobertura de Bosque (%)", "purple")

panel_p <- p4 / p5
ggsave(paste0("~/Desktop/FPVA/Resultados/Aves/Respuesta_Deteccion_", especie_objetivo, ".png"), 
       panel_p, width = 8, height = 10)

# Mostrar en consola
print(panel_psi)
print(panel_p)
