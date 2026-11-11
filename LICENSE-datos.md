# Licencias de los datos y del texto

- **Código** (`src/`, `scripts/`, `pruebas/`, `manuscrito/armar_entregables.py`, `Makefile`): MIT (`LICENSE`). Incluye, sin cambios, el paquete común de Estadísticas Vitales (`src/eevv_insumos/`, `scripts/insumos_eevv/`, `pruebas/insumos_eevv/`) y el clasificador de códigos CIE-10 del repositorio de Estadísticas Vitales (`src/agrupamiento/codigos.py`), ambos de Daniel A. Garavito Jiménez bajo la misma licencia.
- **Texto del manuscrito, tablas y figuras** (`manuscrito/`, `output/`): Creative Commons Atribución-NoComercial 4.0 Internacional (CC BY-NC 4.0), © 2026 Daniel A. Garavito Jiménez; citar a Sinergística con enlace.
- **Documentación de método** (`docs/`, `README.md`, `data/README.md`): Creative Commons Atribución 4.0 Internacional (CC BY 4.0), © 2026 Daniel A. Garavito Jiménez. Los mapas usan la geometría del Marco Geoestadístico Nacional del DANE («Fuente: DANE»).
- **Datos derivados de `data/intermedios/`**: CC BY-NC 4.0, salvo las columnas derivadas de fuentes CC BY-SA, que conservan CC BY-SA 4.0 (abajo). Cada fuente original se cita así:
  - Estadísticas Vitales y proyecciones de población: «Fuente: Departamento Administrativo Nacional de Estadística: www.dane.gov.co». Son agregados municipio × causa (y × año) sin registros individuales, con la misma o menor granularidad que el cuadro 5 de defunciones del DANE; las proporciones municipales de calidad del registro se publican sin los conteos (`docs/INSUMOS.md`, sección 4).
  - CNPV 2018, MGN 2018, NBI, IPM y etnia (DANE, geoportal): uso libre con la cita «Fuente: DANE».
  - DIVIPOLA (DANE), REPS 2022 (Ministerio de Salud y Protección Social) y Catálogo Nacional de Estaciones (IDEAM), de datos.gov.co: CC BY-SA 4.0. Las columnas de `covariables.parquet` derivadas de esas fuentes (coordenadas, altitud, camas, distancias) se comparten con la misma licencia CC BY-SA 4.0.
  - Clasificación PDET: norma pública (Decreto Ley 893 de 2017).
- **Insumos originales de `data/raw/`**: rigen los términos de cada fuente (`data/raw/MANIFIESTO.csv`, columna `licencia`). Los microdatos del DANE no se redistribuyen: los términos del Archivo Nacional de Datos prohíben ponerlos a disposición de múltiples usuarios sin visto bueno escrito del DANE (`docs/INSUMOS_EEVV.md`). El cuadro 5 de defunciones 2018 tampoco se redistribuye; se descarga del portal del DANE.

## Cambio de licencia

El 4 de noviembre de 2026, antes de publicar el repositorio, el texto, las tablas, las figuras y los datos derivados pasaron de CC BY 4.0 a CC BY-NC 4.0, según la política de licencias por capa de Sinergística. La documentación de método sigue en CC BY 4.0 y el código, en MIT.
