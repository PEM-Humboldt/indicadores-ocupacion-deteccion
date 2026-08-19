# ==============================================================================
# 04_exportar_resultados.R
# Grafica la relacion entre psi y las covariables, y exporta los resultados
# a un libro de Excel.
#
# Requiere haber corrido antes los scripts 01 a 03 (en orden).
# ==============================================================================

library(ggplot2)
library(openxlsx)

# --- 1. Curva de respuesta: psi vs. % de bosque (dejando rio y via en su media) --
rango_bosque <- seq(min(covariables_sitio$pct_bosque),
                     max(covariables_sitio$pct_bosque), length.out = 50)

nuevos_datos <- data.frame(
  pct_bosque = estandarizar(rango_bosque, medias["pct_bosque"], sds["pct_bosque"]),
  dist_rio_m = 0,  # covariable en su media (estandarizada = 0)
  dist_via_m = 0
)

pred_curva <- predict(modelo_cov, type = "state", newdata = nuevos_datos)
pred_curva$pct_bosque <- rango_bosque

p_respuesta <- ggplot(pred_curva, aes(x = pct_bosque, y = Predicted)) +
  geom_ribbon(aes(ymin = pmax(0, Predicted - SE), ymax = pmin(1, Predicted + SE)),
              alpha = 0.2, fill = "forestgreen") +
  geom_line(color = "forestgreen", linewidth = 1.1) +
  labs(title = paste("Curva de respuesta —", especie_objetivo),
       subtitle = "Probabilidad de ocupacion (psi) segun % de cobertura boscosa",
       x = "% de bosque en el sitio", y = "psi") +
  ylim(0, 1) +
  theme_bw(base_size = 11)

ggsave(file.path(path_out, "01_curva_respuesta_bosque.png"),
       p_respuesta, width = 8, height = 5, dpi = 150)

# --- 2. Mapa simple del grid de prediccion (dispersion, sin GIS) ------------
p_grid <- ggplot(grid_prediccion, aes(x = dist_via_m, y = pct_bosque,
                                       color = psi_predicho, size = psi_predicho)) +
  geom_point(alpha = 0.85) +
  scale_color_viridis_c(name = "psi", limits = c(0, 1)) +
  labs(title = "Ocupacion predicha por celda del grid",
       subtitle = "Cada punto es una celda del paisaje hipotetico",
       x = "Distancia a vias (m)", y = "% de bosque") +
  theme_bw(base_size = 11)

ggsave(file.path(path_out, "02_grid_ocupacion.png"),
       p_grid, width = 8, height = 5, dpi = 150)

# --- 3. Exportar a Excel -----------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "Indicadores_Ocupacion")
writeData(wb, "Indicadores_Ocupacion", tab_indicadores)

addWorksheet(wb, "Seleccion_Modelo_AIC")
writeData(wb, "Seleccion_Modelo_AIC", tab_aic)

addWorksheet(wb, "Prediccion_Grid")
writeData(wb, "Prediccion_Grid", grid_prediccion)

saveWorkbook(wb, file.path(path_out, "Indicadores_Ocupacion.xlsx"),
             overwrite = TRUE)

cat("=== LISTO ===\n")
cat("Resultados guardados en:", path_out, "\n")
cat("Libro de Excel: Indicadores_Ocupacion.xlsx\n")
cat("Graficos: 01_curva_respuesta_bosque.png, 02_grid_ocupacion.png\n")
