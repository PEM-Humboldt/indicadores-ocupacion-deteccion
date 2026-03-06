# ==========================================================
# 2. PREPARACIÓN DE COVARIABLES (Resampleo a 500m)
# ==========================================================
path_base <- "~/Desktop/FPVA/Analisis/Datos/Geograficos/Co-variables Ocupación/"

# Cargar originales
cov_90m <- c(
  rast(paste0(path_base, "BosqueNoBosque__90mts/agreg_bosque90_2.tif")),
  rast(paste0(path_base, "EucDistance_Drenajes90mts/RiosEudis90.tif")),
  rast(paste0(path_base, "EucliDistance_vias90mts/ViasEudis90.tif"))
)

# --- CAMBIO CLAVE: Agregación a 500m ---
# 500m / 90m ≈ 5.5. Usamos factor 6 para suavizar el grano.
# Para distancias usamos el promedio (mean), para bosque podrías usar el promedio 
# (que se convierte en % de cobertura en ese pixel de 500m).
cov_stack <- aggregate(cov_90m, fact = 6, fun = "mean")
names(cov_stack) <- c("bosque", "dist_rios", "dist_vias")

# Extraer usando los puntos de FOTOTRAMPEO
puntos_sf <- dev_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(cov_stack))

# Extraer valores (ahora representan el promedio en un radio de ~500m)
site_covs_raw <- terra::extract(cov_stack, puntos_sf)
site_covs_raw$deploymentID <- dev_data$deploymentID

# Recalcular parámetros (necesario para las gráficas de respuesta después)
params <- site_covs_raw %>%
  summarise(across(c(bosque, dist_rios, dist_vias), 
                   list(mu = ~mean(.x, na.rm=T), sd = ~sd(.x, na.rm=T))))

site_covs_scaled <- site_covs_raw %>%
  mutate(across(c(bosque, dist_rios, dist_vias), scale)) %>%
  filter(complete.cases(bosque, dist_rios, dist_vias))



###
# Reproyectar

# El resto del código de predicción funcionará igual, pero el raster de salida
# será notablemente más "pixelado" (grano grueso), reflejando los 500m.
env_stack_scaled <- cov_stack # Este ya viene a 500m
env_stack_scaled$bosque    <- (cov_stack$bosque - params$bosque_mu) / params$bosque_sd
env_stack_scaled$dist_rios <- (cov_stack$dist_rios - params$dist_rios_mu) / params$dist_rios_sd
env_stack_scaled$dist_vias <- (cov_stack$dist_vias - params$dist_vias_mu) / params$dist_vias_sd

# ... sigue con el logit_psi ...
#####
# Modelo


plot(env_stack_scaled$dist_vias)
