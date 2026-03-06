
# Indicadores Ocupacion y Deteccion

Este repositorio contiene el flujo de trabajo para el modelamiento de **Ocupación de Sitios (Single-season Occupancy Models)**. El código utiliza modelos jerárquicos para estimar la probabilidad de presencia de especies detectadas mediante fototrampeo y bioacústica, corrigiendo el sesgo por detectabilidad imperfecta y proyectando los resultados espacialmente.

**Estado del Proyecto:** Estable / Finalizado. Incluye la integración de covariables ambientales (Rasters) y la generación de mapas de probabilidad de uso de hábitat.

---

## Prerequisitos

Para ejecutar este análisis, es necesario contar con **R (versión 4.0+)** y las siguientes bibliotecas especializadas en ecología y análisis espacial:

* `unmarked`: Motor estadístico para modelos de ocupación y abundancia.
* `terra` y `sf`: Manejo de objetos ráster y vectores (geoprocesamiento).
* `tidyterra`: Extensión de ggplot2 para visualizar datos espaciales de terra.
* `tidyverse` y `lubridate`: Procesamiento de datos y series de tiempo.

**Instalación:**

```r
install.packages(c("unmarked", "terra", "sf", "tidyterra", "tidyverse", "lubridate", "openxlsx"))

```

---

## Archivos Necesarios

1. **Datos Biológicos:** Archivo Excel con las pestañas `Observations`, `Media` y `Deployment`.
2. **Covariables Ambientales:** Archivos en formato `.tif` (Ráster) que representen variables del paisaje. El script espera:
* `agreg_bosque90_2.tif`: Cobertura boscosa.
* `RiosEudis90.tif`: Distancia euclidiana a drenajes.
* `ViasEudis90.tif`: Distancia euclidiana a vías.



* **Ubicación de Covariables:** Deben estar alojadas en la ruta definida en `path_base`.

---

## Cómo ejecutar

Siga el orden de los bloques para asegurar la construcción correcta del objeto `unmarkedFrame`:

1. **Preparación de Covariables:**
* El script carga el stack de rásters, extrae los valores para cada punto de cámara y **estandariza los datos** (media=0, sd=1) para mejorar la convergencia de los modelos.


2. **Construcción de Historia de Detección:**
* Se genera una matriz de presencia/ausencia agrupada por ocasiones (ej. intervalos de 15 días).
* *Nota:* El script maneja automáticamente los valores `NA` y asegura la coincidencia entre sitios muestreados y covariables.


3. **Ajuste y Selección de Modelos (AIC):**
* Se ejecutan modelos que prueban efectos en la **Ocupación ($\psi$)** y la **Detección ($p$)**.
* Se utiliza el Criterio de Información de Akaike (AIC) para seleccionar el modelo con mejor soporte estadístico.


4. **Generación de Curvas de Respuesta:**
* El script incluye funciones para des-escalar los datos y graficar la respuesta real de la especie ante las variables ambientales (ej. % de bosque o distancia en metros).


5. **Proyección Espacial (Mapeo):**
* Los coeficientes ($\beta$) del modelo ganador se aplican sobre los rásters originales para generar un **Mapa de Probabilidad de Ocupación** continuo.



### Ejemplo de Salida

* **Tabla_AIC.csv:** Ranking de modelos para selección del mejor candidato.
* **Curvas_Respuesta.png:** Gráficas con intervalos de confianza del 95%.
* **Mapa_Ocupacion.tif:** Archivo ráster listo para ser abierto en **QGIS** o **ArcGIS**.

---

## Interpretación de Resultados

* **$\psi$ (Psi):** Probabilidad de que el sitio esté ocupado por la especie.
* **$p$ (Rho):** Probabilidad de detectar a la especie dado que el sitio está ocupado.
* **Círculo Rojo/Blanco en Mapa:** Ubicación de las estaciones de muestreo sobre la probabilidad proyectada.

---

## Autores(as) y contacto

* **Juan C Rey** - *Investigador* - [jrey@humboldt.org.co]

---

## Licencia

Este proyecto está bajo la licencia MIT.

## Agradecimientos

* A la librería **unmarked** por proporcionar las herramientas para el análisis de poblaciones animales.
* Al equipo GIS por la provisión de las capas ambientales estandarizadas.

---