# Indicadores Ocupación y Detección

Este repositorio contiene el flujo de trabajo para el modelamiento de **Ocupación de Sitios (Single-season Occupancy Models)**. El código utiliza modelos jerárquicos para estimar la probabilidad de presencia de especies detectadas mediante fototrampeo y bioacústica, corrigiendo el sesgo por detectabilidad imperfecta y proyectando los resultados espacialmente.

El siguiente flujo de trabajo integra los siguientes indicadores:
- Porcentaje del área de estudio ocupada por la especie: proporción del área donde el modelo predice presencia (ψ ≥ umbral definido).
- Probabilidad de detección (p, sin covariables): qué tan probable es registrar la especie dado que está presente.
- Probabilidad de ocupación (ψ, sin covariables): probabilidad de que un sitio esté ocupado por la especie.
---
## Fundamento metodológico

Esta sección explica el *qué* y el *por qué* de los cálculos. El detalle de calidad e incertidumbre de los indicadores derivados (porcentaje de área ocupada, psi y p sin covariables) vive en la ficha del catálogo de indicadores; aquí va lo necesario para entender y reproducir el código.

- **Enfoque:** modelos jerárquicos de ocupación de una sola estación (MacKenzie et al., 2017), siguiendo además las mejores prácticas de eBird para el manejo de datos de detección imperfecta (Strimas-Mackey et al., 2023). El diseño asume que la presencia/ausencia observada en un sitio resulta de dos procesos: uno biológico (ocupación, ψ) y uno de observación (detección, p), lo que permite corregir el sesgo de detectabilidad imperfecta.
- **Implementación:** paquete `unmarked` de R. Se ajustan modelos nulos (sin covariables) y modelos con covariables ambientales; la selección del modelo más informativo se hace con el Criterio de Información de Akaike (AIC).
- **Cómo se calcula, paso a paso:**
  1. **Curaduría:** consolidación de las tablas de monitoreo comunitario (Observations, Media, Deployment).
  2. **Preparación de covariables:** carga del stack de rásters, extracción de valores por punto de muestreo y estandarización (Z-score).
  3. **Historia de detección:** construcción de la matriz de presencia/ausencia en ocasiones (por ejemplo, de 7 o 15 días, según el método de detección).
  4. **Modelamiento:** ajuste en `unmarked`, tanto del modelo nulo (intercepto) como de modelos con covariables ambientales.
  5. **Extracción de resultados:**
     - Para ψ y p **sin covariables**, se usan los estimados del intercepto del modelo nulo.
     - Para el **porcentaje de área ocupada**, se proyectan los coeficientes β del modelo seleccionado sobre los rásters y se aplica un umbral de ψ (definido en la ficha del indicador).

---

## Estructura del repositorio


```
.
├── data/
│   ├── I2D_FPVA_Fototrampeo_20260219.xlsx      # Observations, Media, Deployment (fototrampeo)
│   └── I2D_FPVA_Bioacustica_20260219.xlsx      # Observations, Media, Deployment (bioacústica)
├── covariables/
│   ├── agreg_bosque90_2.tif                    # Cobertura boscosa
│   ├── RiosEudis90.tif                         # Distancia euclidiana a drenajes
│   └── ViasEudis90.tif                         # Distancia euclidiana a vías
├── scripts/
│   ├── 00_prep_covariables.R                   # Carga y estandarización del stack de rásters
│   ├── Occu_fototrampeo.R                      # Modelos de ocupación para mamíferos (fototrampeo)
│   ├── Occu_bioacustica.R                      # Modelos de ocupación para aves (bioacústica)
│   └── Occu_puntos_conteo.R                    # Modelos de ocupación para aves (puntos de conteo), si aplica
├── outputs/
│   ├── Tabla_AIC.csv
│   ├── Curvas_Respuesta.png
│   └── Mapa_Ocupacion.tif
├── README.md
└── LICENCIA
```

## Descripción de los scripts

| Script | Qué hace | Insumos que usa | Salida principal |
|---|---|---|---|
| `00_prep_covariables.R` | Carga el stack de rásters ambientales, extrae los valores en cada punto de muestreo y estandariza los datos (media=0, sd=1) para mejorar la convergencia de los modelos. | Rásters `.tif` + coordenadas de las estaciones | Tabla de covariables estandarizadas por sitio |
| `Occu_fototrampeo.R` | Ajusta los modelos de ocupación de una sola estación para las especies de mamíferos detectadas por cámaras trampa; selecciona el mejor modelo por AIC y proyecta el mapa de ocupación. | Datos biológicos (fototrampeo) + covariables estandarizadas | `Tabla_AIC.csv`, `Curvas_Respuesta.png`, `Mapa_Ocupacion.tif` |
| `Occu_bioacustica.R` | Igual que el anterior, pero para especies detectadas por bioacústica. Sigue la misma estructura de `Occu_fototrampeo.R`, con los ajustes propios del método (ventana de ocasión, filtros de validación experta). | Datos biológicos (bioacústica) + covariables estandarizadas | `Tabla_AIC.csv`, `Curvas_Respuesta.png`, `Mapa_Ocupacion.tif` (por especie) |
| `Occu_puntos_conteo.R` | Igual estructura, para especies registradas en puntos de conteo. | Datos biológicos (puntos de conteo) + covariables estandarizadas | Igual que los anteriores |

---

## Prerrequisitos

Para ejecutar este análisis es necesario contar con **R (versión 4.0 o superior)** y las siguientes librerías. Se indican las versiones con las que se probó el flujo, para evitar problemas de compatibilidad en ejecuciones futuras:

| Librería | Uso en el proyecto | Versión probada |
|---|---|---|
| `unmarked` | Motor estadístico para modelos de ocupación y abundancia | ‘1.4.1’|
| `terra` | Manejo de objetos ráster (geoprocesamiento) | ‘1.8.60’ |
| `sf` | Manejo de objetos vectoriales (geoprocesamiento) | ‘1.0.21’ |
| `tidyterra` | Extensión de ggplot2 para visualizar datos espaciales de `terra` | ‘0.7.2’ |
| `tidyverse` | Procesamiento de datos | ‘2.0.0’ |
| `lubridate` | Manejo de series de tiempo | ‘1.9.2’ |
| `openxlsx` | Lectura/escritura de Excel | ‘4.2.5.2’ |


**Instalación:**

```r
install.packages(c("unmarked", "terra", "sf", "tidyterra", "tidyverse", "lubridate", "openxlsx"))
```

---


## Archivos necesarios

1. **Datos biológicos:** archivo Excel con las pestañas `Observations`, `Media` y `Deployment`.
2. **Covariables ambientales:** archivos en formato `.tif` (ráster) que representan variables del paisaje. El script usa:
   - `agreg_bosque90_2.tif`: cobertura boscosa.
   - `RiosEudis90.tif`: distancia euclidiana a drenajes.
   - `ViasEudis90.tif`: distancia euclidiana a vías.
3. **Ubicación de covariables:** deben estar alojadas en la ruta definida en `path_base`.

---

## Cómo ejecutar

El flujo está pensado para correrse **de forma secuencial**. Hay **dos puntos donde el usuario debe interactuar**; el resto corre sin intervención una vez definidas las rutas iniciales:

1. **Preparación de covariables** *(requiere acción del usuario: definir `path_base` con la ubicación de los rásters)*
   - El script carga el stack de rásters, extrae los valores para cada punto de cámara y estandariza los datos (media=0, sd=1) para mejorar la convergencia de los modelos.

2. **Construcción de historia de detección** *(se ejecuta sin intervención)*
   - Se genera una matriz de presencia/ausencia agrupada por ocasiones (ej. intervalos de 15 días).
   - *Nota:* el script maneja automáticamente los valores `NA` y asegura la coincidencia entre sitios muestreados y covariables.

3. **Ajuste y selección de modelos (AIC)** *(se ejecuta sin intervención, salvo que quieras cambiar el set de covariables candidatas)*
   - Se ejecutan modelos que prueban efectos en la ocupación (ψ) y la detección (p).
   - Se usa el Criterio de Información de Akaike (AIC) para seleccionar el modelo con mejor soporte estadístico.

4. **Generación de curvas de respuesta** *(se ejecuta sin intervención)*
   - El script incluye funciones para des-escalar los datos y graficar la respuesta real de la especie ante las variables ambientales (ej. % de bosque o distancia en metros).

5. **Proyección espacial (mapeo)** *(requiere acción del usuario: confirmar/ajustar el umbral de ψ usado para el porcentaje de área ocupada)*
   - Los coeficientes (β) del modelo ganador se aplican sobre los rásters originales para generar un mapa de probabilidad de ocupación continuo.

### Ejemplo de salida

* **Tabla_AIC.csv:** ranking de modelos para selección del mejor candidato.
* **Curvas_Respuesta.png:** gráficas con intervalos de confianza del 95%.
* **Mapa_Ocupacion.tif:** archivo ráster listo para ser abierto en QGIS o ArcGIS.

---

## Interpretación de resultados

* **ψ (Psi):** probabilidad de que el sitio esté ocupado por la especie.
* **p:** probabilidad de detectar a la especie dado que el sitio está ocupado.
* **Círculo rojo/blanco en el mapa:** ubicación de las estaciones de muestreo sobre la probabilidad proyectada.

---

## Autores(as) y contacto

* **Juan C Rey** — *Investigador* — [jrey@humboldt.org.co]

## Licencia

Este proyecto está bajo la licencia MIT.

## Agradecimientos

* Al equipo de campo del proyecto **FPVA**  y a la comunidad por la recolección de los datos biológicos.

## Referencias metodológicas

1. MacKenzie, D.I. et al. (2017). *Occupancy Estimation and Modeling: Inferring Patterns and Dynamics of Species Occurrence* (2nd ed.).
2. Royle, J.A., Nichols, J.D. & Kéry, M. (2005). Modelling occurrence and abundance of species when detection is imperfect. *Oikos*, 110: 353–359.
3. Ruiz-Gutiérrez, V., Zipkin, E.F. & Dhondt, A.A. (2010). Occupancy dynamics in a tropical bird community. *Journal of Applied Ecology*, 47: 621–630.
4. Sollmann, R. (2018). A gentle introduction to camera-trap data analysis. *African Journal of Ecology*, 56: 740–749.
5. Strimas-Mackey, M. et al. (2023). Best Practices for Using eBird Data.
