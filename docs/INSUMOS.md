# Insumos del artículo A3: procedencia, términos y validación

Una ficha por cada insumo propio del artículo (grupos «a3» y «a3-config» de `data/raw/MANIFIESTO.csv`). Los microdatos de Estadísticas Vitales, las proyecciones de población y las estimaciones de mortalidad del DANE (grupo «eevv», 79 archivos) se procesan con el paquete común sin cambios y tienen sus fichas en `docs/INSUMOS_EEVV.md`. La huella SHA-256, el tamaño y la fecha de descarga de cada archivo están en el manifiesto; `make verificar-insumos` y `pruebas/test_insumos.py` las comprueban.

## Contenido

1. Resumen de insumos
2. Fichas de los insumos propios
3. Código vendorizado y cálculos tomados de otros trabajos
4. Agregados publicados y regla de celdas pequeñas
5. Validación contra cifras publicadas

## 1. Resumen de insumos

| Insumo | Archivo | Tamaño (B) | ¿Versionado? | Uso |
|---|---|--:|---|---|
| Microdatos de defunciones no fetales 2015–2025 (DANE, ANDA) | `data/raw/*.zip` | | no (términos del ANDA) | `scripts/insumos_eevv/` → `data/processed/ext/def_<año>.parquet`, `defunciones_eevv.parquet` |
| Proyecciones de población municipal y departamental (DANE) | `data/raw/poblacion/` | | no | `scripts/insumos_eevv/04_poblacion.py` → `mun_pob_edad.parquet`, `poblacion_eevv.parquet` |
| Tasas de mortalidad por edad, sexo y departamento 1985–2017 (DANE) | `data/raw/poblacion/dane_mor/DCD-Mor-EstSexNal-Dep-1985-2017_VP.xlsx` | 1.815.803 | no (se descarga) | completitud departamental |
| MGN 2018 integrado con el CNPV 2018, municipios (DANE) | `data/raw/externos/63/mgn2018_mpio_cnpv.geojson` | 3.560.508 | sí | internet, educación superior, área, mapas |
| MGN 2018, municipios, geometría fina (DANE) | `data/raw/externos/63/mgn2018_mpio_fino.geojson` | 8.494.781 | sí | matriz de contigüidad |
| Diccionario del MGN integrado (DANE, 2020) | `data/raw/externos/63/instructivo_mgn_integrado_cnpv2018.pdf` | 591.019 | sí | documentación |
| NBI 2018 por municipio (DANE) | `data/raw/externos/63/nbi2018_mpio.json` | | sí | NBI |
| IPM 2018 por municipio (DANE) | `data/raw/externos/63/ipm2018_mpio.json` | | sí | descriptivo |
| Población indígena y afrocolombiana 2018 (DANE) | `data/raw/externos/63/etnia2018_mpio.json` | | sí | etnia |
| REPS 2022, capacidad instalada, todas las filas (MinSalud) | `data/raw/externos/63/reps_capacidad_2022.csv` | | sí | camas, camas de cuidado intensivo |
| DIVIPOLA con coordenadas (DANE) | `data/raw/externos/52/divipola_gdxc-w37w.csv` | 70.850 | sí | puntos de referencia, distancias, vecinos más cercanos |
| Catálogo Nacional de Estaciones (IDEAM) | `data/raw/externos/52/ideam_estaciones.csv` | 3.455.945 | sí | altitud |
| Cuadro 5 de defunciones 2018 (DANE) | `data/raw/externos/63/CUADRO5-NOFETALES-2018-def.xls` | 8.269.824 | no (se descarga) | validación de los observados y regla de celdas |
| Clasificación PDET | `config/clasificacion_municipal.csv` | 11.293 | sí | PDET |
| Tipología de códigos basura (Naghavi et al., 2010) | `config/codigos_basura_naghavi2010.csv` | | sí | calidad del registro |
| Selección de causas del trabajo previo | `config/seleccion_causas_departamental.csv` | | sí | las 12 causas |
| Cifras citadas de otros trabajos | `config/cifras_citadas.csv` | | sí | verificación de cifras del texto |

Los tamaños y huellas exactos de todos los archivos están en el manifiesto.

## 2. Fichas de los insumos propios

### MGN 2018 integrado con el CNPV 2018

- **Origen:** DANE, Marco Geoestadístico Nacional, versión 2018, nivel municipio integrado con variables anonimizadas del Censo Nacional de Población y Vivienda 2018. Servicio ArcGIS REST `MGN2018/Serv_CapaMunicipiosInt_2018/MapServer/0` del geoportal del DANE (`https://portalgis.dane.gov.co/mparcgis/rest/services/`), consulta con `outFields=*`, `maxAllowableOffset=0.003`, `geometryPrecision=5`, `outSR=4326`, `f=geojson` (URL completa en el manifiesto). Descarga: 24 de octubre de 2026. El servidor no publica `robots.txt`.
- **Archivo:** 1.122 polígonos (1.103 municipios, 18 áreas no municipalizadas y una isla) con 86 atributos; se usan `MPIO_CDPMP`, `AREA`, `LATITUD`, `LONGITUD`, `STP19_INT1/2` (viviendas con y sin internet), `STP51_PRIM`, `STP51_SECU`, `STP51_SUPE`, `STP51_POST`, `STP51_13_E` (nivel educativo) y `STP27_PERS`.
- **Definiciones:** tabla 2 de `instructivo_mgn_integrado_cnpv2018.pdf` (DANE, 2020). Las variables son anonimizadas: la suma de `STP51_*` cubre cerca de la mitad de las personas, por eso la educación superior se usa como proporción entre quienes informan nivel.
- **Términos:** uso libre con la cita «Fuente: DANE» (términos de www.dane.gov.co y del geoportal). Se redistribuye.
- **Transformaciones:** `src/agrupamiento/datos.py:covariables()` (internet = INT1/(INT1+INT2); superior = (SUPE+POST)/(PRIM+SECU+SUPE+POST+NINGUNO); densidad = población media 2018–2025 / área); `figuras.py:figuras_y_tablas()` (mapas).
- **Validación:** 1.122 polígonos y suma de `STP27_PERS` = 44.164.417, la población censada del CNPV 2018 publicada por el DANE (`test_mgn_cnpv_personas_censadas`).

### MGN 2018, geometría fina

- **Origen:** la misma capa, solo `MPIO_CDPMP`, con `maxAllowableOffset=0.0005` y `geometryPrecision=6`. Descarga: 24 de octubre de 2026.
- **Uso:** `src/agrupamiento/pesos.py:reina_desde_mgn()` construye la contigüidad tipo reina (con `shapely.make_valid` y `libpysal.weights.Queen`); la geometría de 0,003 grados deja huecos entre vecinos y por eso se descargó esta. Resultado: 6.396 pares ordenados de vecinos (5,7 por municipio), guardados en `data/intermedios/pesos_reina.csv`.
- **Validación:** los mismos 1.122 códigos que la capa integrada (`test_mgn_cnpv_personas_censadas`); la lista de adyacencia reconstruida en el nivel B es idéntica a la versionada (`scripts/05_comparar_intermedios.py`).

### NBI, IPM y población étnica 2018

- **Origen:** servicios del DANE `INDICADORES_COND_DE_VIDA/Serv_Mpios_NBI_PropPersonasEnNBI_Total_2018/MapServer/4` (campo `NBIC_Total_Prop_de_Personas_en_NBIPorc`), `INDICADORES_COND_DE_VIDA/Serv_Mpios_IndPobrezaMultidimensional_2018/MapServer/4` (`IPM`) e `INDICADORES_GRUPOS_ETNICOS/Serv_Mpios_PoblacionIndigena_2018/MapServer/4` (personas y porcentaje indígena y negro, mulato, afrodescendiente o afrocolombiano), enlazados desde datos.gov.co. Las capas de NBI y etnia se bajaron en dos páginas de 1.000 registros y se unieron. Descarga: 24 de octubre de 2026.
- **Términos:** uso libre con la cita «Fuente: DANE». Se redistribuyen.
- **Transformaciones:** `datos.py:covariables()`. Las celdas vacías de etnia corresponden a municipios sin personas del grupo (también viene vacío el total de unidades) y se leen como cero: 152 municipios en la indígena y 30 en la afrocolombiana.
- **Validación:** 1.122 municipios en cada capa; NBI e IPM entre 0 y 100 (`test_nbi_ipm`); suma de la población indígena = 1.905.617 y de la afrocolombiana = 2.950.072, las cifras del CNPV 2018 que publica el DANE en sus presentaciones de grupos étnicos de 2019 (`test_etnia_totales_publicados`).

### REPS 2022, capacidad instalada

- **Origen:** Ministerio de Salud y Protección Social, «Relación de IPS públicas y privadas según el nivel de atención y capacidad instalada», datos.gov.co `s2ru-bqt6`, siete columnas, todas las filas. Corte REPS: 5 de noviembre de 2022. Descarga: 24 de octubre de 2026.
- **Archivo:** 41.427 filas.
- **Términos:** CC BY-SA 4.0 (`meta_s2ru.json`). Se redistribuye.
- **Transformaciones:** `datos.py:covariables()`: el municipio sale de las cinco primeras cifras del código de sede completado a diez dígitos; se suman las camas y las camas cuya descripción contiene «intensiv» (93 municipios las tienen).
- **Validación:** 41.427 filas, igual al conteo que informa la API (`count(*)`), y la fecha de corte en todas (`test_reps_2022`).

### DIVIPOLA con coordenadas y catálogo de estaciones del IDEAM

- **Origen:** los mismos archivos del estudio 52 de la serie (fichas completas en el repositorio del artículo A2): DANE, DIVIPOLA, datos.gov.co `gdxc-w37w`, corte del 30 de diciembre de 2024; IDEAM, Catálogo Nacional de Estaciones, datos.gov.co `hp9r-jxuu`. Descarga: 23 de octubre de 2026 (el catálogo del IDEAM se actualizó ese día; una descarga nueva puede tener otra huella).
- **Términos:** CC BY-SA 4.0. Se redistribuyen.
- **Transformaciones:** `datos.py:covariables()`: punto de referencia del municipio (para Mapiripana, 94663, ausente de la DIVIPOLA de 2024, el centroide del MGN 2018); altitud de la estación del IDEAM más cercana (mediana de la distancia: 0,8 km); distancias en línea recta (fórmula del haversine) a la capital del departamento y al municipio más cercano con camas de cuidado intensivo.
- **Validación:** 1.122 entidades y 9.914 estaciones (`test_divipola_ideam_pdet`).

### Cuadro 5 de defunciones 2018

- **Origen:** DANE, Estadísticas Vitales, cifras definitivas 2018, cuadro 5: «Defunciones por grupos de edad y sexo, según departamento, municipio de residencia y grupos de causas de defunción (lista de causas agrupadas 6/67 CIE-10 de OPS)». URL en el manifiesto. Descarga: 29 de octubre de 2026.
- **Términos:** publicación del portal del DANE con la cita «Fuente: DANE». No se redistribuye; lo baja `scripts/00_descargar_insumos.py`.
- **Uso:** solo en `pruebas/test_insumos.py:test_observados_2018_contra_cuadro5`, que valida los observados y sostiene la regla de celdas (sección 4).

### Clasificación PDET

- **Origen:** Agencia de Renovación del Territorio, datos.gov.co `idrk-ba8y` (Decreto Ley 893 de 2017); compilado en `config/clasificacion_municipal.csv` con su nota `clasificacion_municipal.FUENTES.md`. Descarga: 3 de octubre de 2026.
- **Validación:** 170 municipios, como el decreto (`test_divipola_ideam_pdet`).

### Tipología de códigos basura

- **Origen:** Naghavi et al. (2010), *Population Health Metrics* 8:9, cuadro 2, transcrita en el estudio 55 de la serie. Las notas de transcripción están en la columna `nota`.
- **Uso:** `datos.py:calidad_municipal()` y `basura_sin_hipertension()`, con el clasificador de `src/agrupamiento/codigos.py`.

### Selección de causas del trabajo previo

- **Origen:** Garavito Jiménez (2026a), «Dónde el riesgo se agrupa». Las 61 causas con 500 o más defunciones en 2020–2024, su I de Moran departamental de 2022, el valor p con 9.999 permutaciones y el q de Benjamini-Hochberg, calculados con `scripts/10_espacial_667.py` y `figura_p3_causas_v2.py` del repositorio de Estadísticas Vitales sobre `output/tablas/espacial_667.json` (SHA-256 `86cea80c21b325f1c71dec538e776571d5f2c93db3fb9663aa4d899b7e22f5ae`).
- **Uso:** define las 12 causas del artículo (pretest). Es un insumo fijo: su recálculo exige el repositorio de Estadísticas Vitales. Este repositorio sí recalcula la I departamental de 2022 de las 14 series con su propia armonización (`scripts/11_departamental.py`), y encuentra que tres de las 12 causas no se agrupan entre departamentos (`docs/REVISION.md`, sección 2).
- **Validación:** 61 causas, 27 con p < 0,05 y 12 con q < 0,05 (`test_seleccion_causas`).

### Cifras citadas

`config/cifras_citadas.csv` lista las cinco cifras que el manuscrito toma del estudio 55 de la serie (certificación, asistencia médica, códigos basura de 2023–2025 y la brecha de completitud), con el lugar exacto de la fuente. No se recalculan aquí; `verificar_cifras.py` comprueba que el texto las cite igual.

## 3. Código vendorizado y cálculos tomados de otros trabajos

- **Paquete común de Estadísticas Vitales** (`src/eevv_insumos/`, `scripts/insumos_eevv/`, `pruebas/insumos_eevv/`): copia sin cambios de `articulos-cientificos/_comun/`. Ver `docs/INSUMOS_EEVV.md`, `docs/DICCIONARIO_EEVV.md` y `docs/REPLICACION_EEVV.md`.
- **Indicadores municipales de calidad del registro (estudio 55):** `datos.py:calidad_municipal()` reimplementa, con la base extendida del paquete común, las definiciones de `scripts/55_calidad_dato_defunciones.py` del repositorio de Estadísticas Vitales (muertes por causas no externas sin asistencia, códigos basura de los tipos 1 a 4 sobre las muertes con causa CIE-10, muertes en casa). El resultado es idéntico al de ese estudio en los 1.120 municipios (comparación del 30 de octubre de 2026, `docs/REPLICACION.md`).
- **Completitud departamental:** `datos.py:completitud_departamental()` reimplementa la regla de `scripts/35_completitud_cnpv2018.py` del mismo repositorio (registradas con sexo conocido / esperadas con las tasas del DANE por edad simple, sexo y departamento, 2015–2017), de donde la tomaba el estudio 55. Idéntica en los 32 departamentos de la muestra.
- **Clasificador de códigos CIE-10:** `src/agrupamiento/codigos.py`, copia sin cambios de `src/eevv/codigos.py` del repositorio de Estadísticas Vitales.

## 4. Agregados publicados y regla de celdas pequeñas

El estándar pide suprimir las celdas menores de 5 en cruces sensibles. Los agregados del nivel A tienen celdas pequeñas: de las 15.680 celdas municipio × causa de `observados_esperados.parquet`, 4.591 tienen entre 1 y 4 defunciones, y de las 125.440 de `observados_esperados_anual.parquet`, 44.048.

**Observados por municipio y causa: no se suprimen.** El DANE publica cruces más finos. Su cuadro 5 de defunciones da, para cada municipio de residencia y cada causa de la lista 6/67, las muertes por grupo de edad y sexo de cada año, con celdas de 1, 2 y 3 casos. Los observados de este repositorio son sumas de esas celdas (sin edad ni sexo; anuales en el panel y de ocho años en lo demás). La prueba `test_observados_2018_contra_cuadro5` comprueba que los de 2018 coinciden con el cuadro en los 1.120 municipios y las 14 series, salvo una celda: Armenia, todas las causas, 2.222 frente a 2.223, por una defunción sin edad que el artículo excluye. Los esperados son productos de tasas nacionales por la población proyectada, que el DANE también publica. Suprimir los observados haría imposible el suavizado, el modelo de conteo y todas las cifras del texto.

**Covariables de calidad del registro: se publican como proporciones, sin los conteos.** Las muertes sin asistencia médica, con código basura o en casa por municipio cruzan el municipio con una característica del certificado que el DANE publica solo por departamento (cuadro 6, sitio de defunción por departamento de ocurrencia). En los ocho años, el numerador está entre 1 y 4 en 15 municipios (sin asistencia), 20 (códigos basura) y 7 (en casa). Son proporciones sobre todas las muertes del municipio en ocho años, sin edad, sexo, causa ni pertenencia étnica; no describen personas y su denominador es el total de defunciones, que el DANE publica por municipio. Por eso se publican sin supresión. **Esta decisión queda para revisión de Daniel:** si se exige suprimir, la alternativa publicable es reemplazar esos 42 valores por la mediana del departamento; eso cambiaría las estimaciones y el nivel A ya no reproduciría las cifras exactas del manuscrito (habría que recalcularlas con los valores suprimidos y declarar la diferencia).

**Lo que no se publica.** Ningún archivo con una fila por defunción: las bases del paquete común (`data/processed/`) y los conteos por edad (`data/interim/`) no se versionan.

## 5. Validación contra cifras publicadas

| Prueba | Qué compara | Resultado |
|---|---|---|
| `test_intermedios_y_defunciones_publicadas` | defunciones de cada año 2018–2025 frente al total del boletín técnico del DANE; Σ O = Σ E por causa | iguales en los ocho años; 2.231.854 analizadas |
| `test_observados_2018_contra_cuadro5` | observados de 2018 por municipio y causa frente al cuadro 5 del DANE | iguales salvo una celda (una defunción sin edad) |
| `test_mgn_cnpv_personas_censadas` | personas del MGN integrado frente a la población censada del CNPV 2018 | 44.164.417 en ambos |
| `test_etnia_totales_publicados` | población indígena y afrocolombiana frente a las presentaciones del DANE | 1.905.617 y 2.950.072 en ambos |
| `test_reps_2022` | filas y fecha de corte del REPS | 41.427 filas, 5 de noviembre de 2022 |
| `test_divipola_ideam_pdet` | entidades, estaciones y municipios PDET | 1.122, 9.914 y 170 |
| `test_huellas_versionados` | huella y tamaño de cada insumo versionado | todas coinciden |
| `pruebas/insumos_eevv/` | huellas de los zips del DANE, columnas y totales anuales del paquete común | ver `docs/REPLICACION_EEVV.md` |
