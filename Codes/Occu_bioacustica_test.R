library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(tidyterra)
library(openxlsx) 
library(data.table)
# ==========================================================
# 1. CONFIGURACIÓN DE DATOS Y ESPECIE
# ==========================================================

# Cambiamos las referencias a los objetos de fototrampeo
path_data <- "~/Desktop/FPVA/Data/Bioacustica/Plantilla monitoreo acústico - FPV Amazonía.xlsx"
path_out  <- "~/Desktop/FPVA/Resultados/Bioacustica//"
if(!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

obs_data<- fread("~/Desktop/FPVA/Data/Bioacustica/observations.csv")
dev_data  <- fread("~/Desktop/FPVA/Data/Bioacustica/deployments.csv")
med_data <- read.xlsx(path_data, sheet = "Media")
#dev_data  <- read.xlsx(path_data, sheet = "Deployment")

especie_objetivo <- "Lipaugus vociferans"

# Cambiamos las referencias a los objetos de fototrampeo
#obs_data <- bioacustic_obs
#dev_data <- bioacustic_dev
#med_data <- bioacustic_media
# ==========================================================
# 2. PREPARACIÓN DE COVARIABLES
# ==========================================================

# --- Rutas ---
path_bos  <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/Co-variables ocupación/BosqueNoBosque__90mts"
path_dens <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/VariablesOcupacion_Densidad/"
path_dist <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/VariablesOcupacion_2/"
path_dist <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/"

resolucion_m <- 100

# --- Cargar capas ---
r_bosque    <- rast(file.path(path_bos,  "agreg_bosque90_2.tif"))
r_dens_dren <- rast(file.path(path_dens, "dens_drenajes_100m.tif"))
r_dens_vias <- rast(file.path(path_dens, "dens_vias_100m.tif"))
r_dist_dren <- rast(file.path(path_dist, "DrenajesCompleto_2.tif/DrenajesCompleto_2.tif"))
r_dist_vias <- rast(file.path(path_dist, "viascompletas_2/viascompletas_2.tif"))
r_bosque_VA <- rast(file.path(path_dist, "viascompletas_2/viascompletas_2.tif"))



cat("Bosque cargado. Rango:",
    round(minmax(r_bosque)[1], 2), "a",
    round(minmax(r_bosque)[2], 2), "\n")
cat("CRS bosque:", crs(r_bosque, describe = TRUE)$name, "\n")

# --- Función de alineación ---
alinear_capa <- function(r, referencia) {
  if (!same.crs(r, referencia)) r <- project(r, crs(referencia))
  r <- crop(r, ext(referencia))
  r <- resample(r, referencia, method = "bilinear")
  return(r)
}

# --- 2.1 Usar bosque como referencia de grid ---
r_ref <- r_bosque
names(r_ref) <- "pct_bosque"

# --- 2.2 Alinear densidades al grid del bosque ---
r_dens_dren_a        <- alinear_capa(r_dens_dren, r_ref)
r_dens_vias_a        <- alinear_capa(r_dens_vias, r_ref)
r_dist_dren_a        <- alinear_capa(r_dist_dren, r_ref)
r_dist_vias_a        <- alinear_capa(r_dist_vias, r_ref)
r_bosque_VA_a        <- alinear_capa(r_bosque_VA, r_ref)

names(r_dens_dren_a) <- "dens_drenajes"
names(r_dens_vias_a) <- "dens_vias"
names(r_dens_dren_a) <- "dist_drenajes"
names(r_dens_vias_a) <- "dist_vias"
names(r_dens_dren_a) <- "Bosque_VA"

# --- 2.3 Stack final: 3 covariables ---
cov_stack <- c(r_ref, r_dens_dren_a, r_dens_vias_a, r_dist_dren_a, r_dist_vias_a, r_bosque_VA_a)
cov_names <- c("pct_bosque", "dens_drenajes", "dens_vias", "dist_drenajes", "dist_vias", "Bosque_VA")
names(cov_stack) <- cov_names

cat("¿Stack alineado?",
    compareGeom(cov_stack[[1]], cov_stack[[3]], stopOnError = FALSE), "\n")
cat("Covariables:", paste(cov_names, collapse = ", "), "\n")

# --- 2.4 Extracción en puntos de muestreo ---
puntos_sf <- dev_data |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(cov_stack))

site_covs_raw <- terra::extract(cov_stack, puntos_sf)
site_covs_raw$deploymentID <- dev_data$deploymentID

params <- site_covs_raw |>
  summarise(across(all_of(cov_names),
                   list(mu = ~mean(.x, na.rm = TRUE),
                        sd  = ~sd(.x,   na.rm = TRUE))))

site_covs_scaled <- site_covs_raw |>
  mutate(across(all_of(cov_names), scale)) |>
  filter(complete.cases(pick(all_of(cov_names))))


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

site_covs_unmarked <- site_covs_scaled %>% select(pct_bosque, dens_drenajes, dens_vias, dist_drenajes, dist_vias, Bosque_VA)
umf <- unmarkedFrameOccu(y = y_final, siteCovs = site_covs_unmarked)

# Ajuste de modelos
fm0 <- occu(~1 ~1, data = umf)
fm1 <- occu(~1 ~NDVI, data = umf)
fm2 <- occu(~1 ~NDWI, data = umf)
fm3 <- occu(~1 ~dist_rios, data = umf)
fm4 <- occu(~1 ~dist_rios_1, data = umf)
fm5 <- occu(~1 ~dist_rios_2, data = umf)
fm6 <- occu(~1 ~dist_vias, data = umf)
fm7 <- occu(~1 ~dist_vias_s, data = umf)
fm8 <- occu(~1 ~dist_vias_p, data = umf)
fm9 <- occu(~1 ~NDVI + NDWI, data = umf)
fm10 <- occu(~1 ~dist_r_casas, data = umf)
fm_global <- occu(~1 ~ NDVI + NDWI + dist_rios + dist_rios_1 + dist_rios_2 + dist_vias + 
                    dist_vias_s + dist_vias_p , data = umf)

model_list <- fitList('Null'=fm0, 'NDVI'=fm1, 'NDWI'=fm2, 'Rios'=fm3, 
                      'Rios_1'=fm4, 'Rios_2'=fm5, 'Vias'=fm6,
                      'Vias_sec'=fm7, 'Vias_p'=fm8,'AYV'=fm9, 'Casa'=fm10,  'Global'=fm_global)
print(modSel(model_list))

# Extraer la tabla de selección como data.frame correctamente
tabla_aic <- as(modSel(model_list), "data.frame")

# Exportar ahora sí sin errores
write.csv(tabla_aic, 
          "~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Nothocrax_urumutum.csv", 
          row.names = FALSE)


# ==========================================================
# 4. SELECCION DEL MODELO
# ==========================================================


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

# Guardar curva de NDVI
p_NDVI <- plot_occu_tapir(fm8, "NDVI", "NDVI (%)", "darkgreen")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_NDVI_Nothocrax_urumutum.png", 
       plot = p_NDVI, width = 7, height = 5)

# Guardar curva de NDWI
p_NDWI <- plot_occu_tapir(fm2, "NDWI", "NDWI", "royalblue")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_NDWI_Nothocrax_urumutum.png", 
       plot = p_NDWI, width = 7, height = 5)

# Guardar curva de Ríos
p_rios <- plot_occu_tapir(fm3, "dist_rios", "Distancia a Ríos (m)", "blue")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_rios_Nothocrax_urumutum.png", 
       plot = p_rios, width = 7, height = 5)

# Guardar curva de Ríos
p_rios_s <- plot_occu_tapir(fm4, "dist_rios_1", "Distancia a drenajes simples (m)", "skyblue")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_DrenajesDobles_Nothocrax_urumutum.png", 
       plot = p_rios, width = 7, height = 5)

# Guardar curva de Ríos
p_rios_d <- plot_occu_tapir(fm5, "dist_rios_2", "Distancia a drenajes dobles (m)", "darkblue")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_dreanajesSimple_Nothocrax_urumutum.png", 
       plot = p_rios, width = 7, height = 5)


# Guardar curva de Vías
p_vias <- plot_occu_tapir(fm6, "dist_vias", "Distancia a Vías (m)", "darkred")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_Vias_Nothocrax_urumutum.png", 
       plot = p_vias, width = 7, height = 5)

# Guardar curva de Vías primarias
p_vias_s <- plot_occu_tapir(fm7, "dist_vias_s", "Distancia a Vías Secundarias (m)", "firebrick")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_Vias_secundarias_Nothocrax_urumutum.png", 
       plot = p_vias, width = 7, height = 5)

# Guardar curva de Vías secundarias
p_vias_p <- plot_occu_tapir(fm8, "dist_vias_p", "Distancia a Vías Primarias (m)", "tomato")
ggsave("~/Desktop/FPVA/Resultados/Bioacustica/Nothocrax_urumutum/Curva_Vias_primarias_Nothocrax_urumutum.png", 
       plot = p_vias, width = 7, height = 5)


# ==========================================================
# 5. MAPA DE OCUPACIÓN (Proyección espacial)
# ==========================================================
# Seleccionamos el mejor modelo (supongamos el global)
# Incluir el model avarage para proyectarlo

modelo_ganador <- fm6 
beta <- coef(modelo_ganador, type = "state")

# Escalar el stack para predicción
env_stack_scaled <- cov_stack
env_stack_scaled$NDVI    <- (cov_stack$NDVI - params$NDVI_mu) / params$NDVI_sd
env_stack_scaled$NDWI <- (cov_stack$NDWI - params$NDWI_mu) / params$NDWI_sd
env_stack_scaled$dist_rios <- (cov_stack$dist_rios - params$dist_rios_mu) / params$dist_rios_sd
env_stack_scaled$dist_rios_1    <- (cov_stack$dist_rios_1 - params$dist_rios_1_mu) / params$dist_rios_1_sd
env_stack_scaled$dist_rios_2 <- (cov_stack$dist_rios_2 - params$dist_rios_2_mu) / params$dist_rios_2_sd
env_stack_scaled$dist_vias <- (cov_stack$dist_vias - params$dist_vias_mu) / params$dist_vias_sd
env_stack_scaled$dist_vias_s <- (cov_stack$dist_vias_s - params$dist_vias_s_mu) / params$dist_vias_s_sd
env_stack_scaled$dist_vias_p <- (cov_stack$dist_vias_p - params$dist_vias_p_mu) / params$dist_vias_p_sd

# Cálculo de probabilidad
logit_psi <- beta[1] + 
  beta[2] * env_stack_scaled$NDVI + 
  beta[3] * env_stack_scaled$NDWI + 
  beta[4] * env_stack_scaled$dist_rios +
  beta[5] * env_stack_scaled$dist_rios_1 + 
  beta[6] * env_stack_scaled$dist_rios_2 + 
  beta[7] * env_stack_scaled$dist_vias +
  beta[8] * env_stack_scaled$dist_vias_p +
  beta[9] * env_stack_scaled$dist_vias_p

  

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
