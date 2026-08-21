# Indicadores de ocupación y detección

Este repositorio calcula **3 indicadores de ocupación** a partir de datos de
monitoreo con visitas repetidas (fototrampeo, bioacústica, puntos de conteo,
etc.), usando modelos de ocupación de una sola temporada (MacKenzie et al.
2017). Está pensado para que cualquier persona con conocimientos básicos de R
pueda ejecutarlo de principio a fin con **solo tres covariables ambientales
sencillas**, sin necesidad de manejar archivos ráster ni herramientas de SIG.

## Indicadores que calcula

| # | Indicador | Qué significa |
|---|---|---|
| 1 | Probabilidad de ocupación (ψ, sin covariables) | Probabilidad de que un sitio esté ocupado por la especie |
| 2 | Probabilidad de detección (p, sin covariables) | Probabilidad de detectar la especie dado que el sitio está ocupado |
| 3 | Porcentaje del área ocupada | % de un paisaje hipotético (grid) donde el modelo predice ψ ≥ 0.5 |

Los tres se calculan con `scripts/03_modelo_ocupacion.R`.

## Simplificación frente a la versión anterior

La versión anterior de este flujo trabajaba con archivos ráster (`.tif`),
reproyecciones y alineación de grids a 100 m, lo que exigía conocimientos de
SIG. Esta versión reemplaza los ráster por una **tabla de covariables por
sitio** y una **tabla de "grid" hipotético** (un CSV con muchas celdas
imaginarias del paisaje) — así el indicador 3 (% de área ocupada) se calcula
igual, pero sin depender de `terra`, `sf` ni de archivos geográficos.

Si más adelante quieres volver a usar covariables ráster reales, basta con
reemplazar `data/grid_prediccion.csv` por los valores extraídos de tus
ráster (una fila por celda) — el resto del flujo no cambia.

## Estructura del repositorio

```
.
├── data/
│   ├── observaciones.csv          # historia de detección por visita (datos de ejemplo)
│   ├── sitios_covariables.csv     # las 3 covariables por sitio de muestreo
│   └── grid_prediccion.csv        # paisaje hipotético para proyectar el % de área ocupada
├── scripts/
│   ├── 00_correr_todo.R           # corre todo el flujo en orden
│   ├── 01_cargar_datos.R
│   ├── 02_historia_deteccion.R
│   ├── 03_modelo_ocupacion.R
│   └── 04_exportar_resultados.R
├── resultados/                     # aquí se guardan el Excel y los gráficos
└── README.md
```

## Datos de ejemplo (inventados)

El repositorio incluye datos simulados de 12 cámaras trampa con 5 visitas cada
una, para una especie de ejemplo (*Dicotyles tajacu*).

**`data/observaciones.csv`** — una fila por visita:

| Columna | Descripción |
|---|---|
| `sitio` | Identificador del sitio (ej. `C01`) |
| `visita` | Número de la ocasión de muestreo (1 a 5) |
| `especie` | Nombre científico |
| `detectado` | `1` si hubo detección en esa visita, `0` si no |

**`data/sitios_covariables.csv`** — las 3 covariables por sitio:

| Columna | Descripción |
|---|---|
| `sitio` | Identificador del sitio |
| `pct_bosque` | % de cobertura boscosa alrededor del sitio |
| `dist_rio_m` | Distancia al río o drenaje más cercano (metros) |
| `dist_via_m` | Distancia a la vía más cercana (metros) |

**`data/grid_prediccion.csv`** — paisaje hipotético para proyectar el modelo:

| Columna | Descripción |
|---|---|
| `id_celda` | Identificador de la celda del paisaje |
| `pct_bosque`, `dist_rio_m`, `dist_via_m` | Mismas 3 covariables, para muchas celdas imaginarias |

Para usar tus propios datos, reemplaza estos tres archivos manteniendo las
mismas columnas. Puedes cambiar la especie a modelar editando
`especie_objetivo` en `scripts/02_historia_deteccion.R`.

## Prerrequisitos

R (4.0 o superior) y las siguientes librerías:

```r
install.packages(c("tidyverse", "unmarked", "openxlsx"))
```

| Librería | Uso en el proyecto | Versión probada |
|---|---|---|
| `unmarked` | Motor estadístico para modelos de ocupación y abundancia | ‘1.4.1’|
| `terra` | Manejo de objetos ráster (geoprocesamiento) | ‘1.8.60’ |
| `sf` | Manejo de objetos vectoriales (geoprocesamiento) | ‘1.0.21’ |
| `tidyterra` | Extensión de ggplot2 para visualizar datos espaciales de `terra` | ‘0.7.2’ |
| `tidyverse` | Procesamiento de datos | ‘2.0.0’ |
| `lubridate` | Manejo de series de tiempo | ‘1.9.2’ |
| `openxlsx` | Lectura/escritura de Excel | ‘4.2.5.2’ |

## Cómo ejecutar

1. Descarga o clona el repositorio y ábrelo como tu directorio de trabajo en R.
2. Corre el script maestro:

   ```r
   source("scripts/00_correr_todo.R")
   ```

3. Revisa la carpeta `resultados/`: encontrarás el libro
   `Indicadores_Ocupacion.xlsx` (con los 3 indicadores, la tabla de selección
   de modelo por AIC y la predicción por celda del grid) y dos gráficos
   (`01_curva_respuesta_bosque.png`, `02_grid_ocupacion.png`).

Si prefieres ejecutar paso a paso, corre los scripts de la carpeta `scripts/`
en el orden en que están numerados (01 a 04); cada uno depende de los
objetos creados por el anterior.

### Cosas que puedes ajustar

- **`especie_objetivo`** en `02_historia_deteccion.R`: la especie a modelar.
- **`UMBRAL_PSI`** en `03_modelo_ocupacion.R`: el umbral de ψ para calcular el
  % de área ocupada (por defecto 0.5).
- El modelo con covariables usa las 3 juntas (`psi ~ pct_bosque + dist_rio_m +
  dist_via_m`); si quieres comparar más combinaciones de covariables, agrega
  más modelos `occu()` a la tabla de selección por AIC en
  `03_modelo_ocupacion.R`.

## Interpretación de resultados

* **ψ (Psi):** probabilidad de que el sitio esté ocupado por la especie.
* **p:** probabilidad de detectar la especie dado que el sitio está ocupado.
* **% de área ocupada:** proporción de celdas del grid de predicción donde
  ψ ≥ el umbral definido.

## Autores(as) y contacto

* **Juan C Rey** — *Investigador* — jrey@humboldt.org.co

## Licencia

Este proyecto está bajo la licencia MIT.

## Referencias metodológicas

1. MacKenzie, D.I. et al. (2017). *Occupancy Estimation and Modeling:
   Inferring Patterns and Dynamics of Species Occurrence* (2nd ed.).
2. Royle, J.A., Nichols, J.D. & Kéry, M. (2005). Modelling occurrence and
   abundance of species when detection is imperfect. *Oikos*, 110: 353–359.
3. Sollmann, R. (2018). A gentle introduction to camera-trap data analysis.
   *African Journal of Ecology*, 56: 740–749.
