# Diccionario de los agregados y de las salidas

Variables de `data/intermedios/` (lo que alimenta el nivel A) y de `output/`. Las bases del paquete común están en `docs/DICCIONARIO_EEVV.md`.

## Contenido

1. data/intermedios/
2. output/
3. Códigos de causa

## 1. data/intermedios/

### observados_esperados.parquet (15.680 filas: 1.120 municipios × 14 series)

| Variable | Tipo | Descripción |
|---|---|---|
| `mun` | entero | código DIVIPOLA del municipio de residencia (Nuevo Belén de Bajirá, 27493, sumado a Riosucio, 27615) |
| `causa` | texto | `todas` o código de la lista 6/67 (sección 3) |
| `O` | entero | defunciones no fetales 2018–2025 con municipio y edad conocidos |
| `E` | real | defunciones esperadas por estandarización indirecta: Σ_a λ_a P_a, con λ_a la tasa de la muestra en el grupo de edad a (0, 1–4, 5–9, …, 80–84, 85 y más) y P_a la población acumulada 2018–2025 |

`observados_esperados_sin_pandemia.parquet`: igual, sin 2020 y 2021. `observados_esperados_anual.parquet` (125.440 filas): igual por año (`anio`), con las tasas de cada año.

### covariables.parquet (1.120 filas)

| Variable | Descripción | Fuente |
|---|---|---|
| `mun`, `dp`, `nmun`, `ndp`, `macro` | municipio, departamento, nombres, macrorregión (Andina, Caribe, Pacífica, Orinoquía, Amazonía, por departamento) | DIVIPOLA |
| `lat`, `lon`, `coord_mgn` | punto de referencia del municipio (grados); `coord_mgn` = el punto es el centroide del MGN porque falta en la DIVIPOLA | DIVIPOLA, MGN |
| `altitud`, `km_estacion` | altitud (m) de la estación del IDEAM más cercana y distancia a ella (km) | IDEAM |
| `area_km2`, `lat_c`, `lon_c`, `personas_cnpv` | área, centroide y personas censadas | MGN 2018 integrado |
| `internet` | viviendas con internet / (con + sin), % | CNPV 2018 |
| `superior` | personas con educación técnica, tecnológica, universitaria o posgrado / personas con nivel informado, % | CNPV 2018 |
| `nbi`, `ipm` | personas con NBI, %; índice de pobreza multidimensional, % | CNPV 2018 |
| `indigena`, `afro` | población indígena y negra, mulata, afrodescendiente o afrocolombiana, % | CNPV 2018 |
| `camas`, `camas_uci` | camas y camas de cuidado intensivo habilitadas en sedes del municipio | REPS 2022 |
| `pdet` | 1 si el municipio es PDET | Decreto Ley 893 de 2017 |
| `def_total` | defunciones 2018–2025 por residencia (todas, incluso sin edad) | EEVV |
| `sin_asistencia` | muertes por causas no externas sin asistencia médica / muertes por causas no externas, %, 2018–2025 | EEVV (definición del estudio 55) |
| `basura` | muertes con código basura (tipos 1–4, Naghavi et al., 2010) / muertes con causa CIE-10, % | EEVV |
| `basura_sin_hta` | ídem entre las muertes con causa distinta de I10–I15 | EEVV |
| `casa` | muertes en casa o domicilio / todas, % | EEVV |
| `completitud` | completitud del registro de defunciones del departamento, 2015–2017, % (registradas / esperadas con las tasas del DANE) | EEVV y DANE |
| `rural` | población en centros poblados y rural disperso, %, 2018–2025 | proyecciones DANE |
| `pob_media` | población media anual 2018–2025 | proyecciones DANE |
| `dist_capital`, `dist_uci` | distancia en línea recta (km) a la capital del departamento y al municipio más cercano con camas de cuidado intensivo | DIVIPOLA, REPS |
| `ldens`, `ldist_capital`, `ldist_uci`, `lcamas` | log(pob_media/área), log(1 + distancias), log(1 + camas por 1.000 habitantes) | derivadas |

### Otros

- `departamental_2022.parquet` (448 filas): `dp`, `causa`, `O`, `E` de 2022 por departamento de residencia (32, sin San Andrés), con la misma estandarización.
- `pesos_reina.csv`: pares (`mun`, `vecino`) de la contigüidad tipo reina sobre el MGN 2018 (6.396 pares ordenados). `pesos_departamentos.csv`: pares (`dp`, `vecino`) derivados.
- `datos_meta.json`: por año, defunciones totales, sin municipio, sin edad (con municipio) y de San Andrés y Providencia; número de municipios.

## 2. output/

| Archivo | Contenido | Etapa |
|---|---|---|
| `resultados.json` | por causa: `moran_y`, `moran_ar` (Assunção-Reis), `I_resid` (por subconjunto de bloques y matriz), `shapley`, `ols` (con LM), `ols_fe`, `conley`, `modelos` (SAR, SEM, SDM, SDEM, SLX), `lr`, `impactos`, `gmm_sdem`, `otras_W`, `regimenes`, `sin_pandemia`, `fe_espacial` | 06, 07 |
| `revision.json` | `EI_H0`; por causa: caída real y placebo, exceso, valor p, numerador y denominador, Shapley con placebo por bloque, pasos no monótonos; `slx` (MCO, MCP, MCP sin I10–I15, con errores agrupados y de Conley, p y q) | 08 |
| `bym2.json` | por causa: modelo nulo y completo (σ², φ, frontera, β, covarianza de β, verosimilitud marginal), I del efecto aleatorio; para 302, `sin_hta` | 09 |
| `panel.json` | por causa: `estimable`, γ, δ y errores agrupados con efectos de año o de departamento × año (`null` si no es estimable: la leucemia) | 10 |
| `departamental.json` | por causa: I de Moran departamental de 2022, p con 999 permutaciones, I del trabajo previo | 11 |
| `tablas/tabla1_covariables.csv` … `tabla10_panel.csv` | tablas del manuscrito y del anexo (tabla 7: placebo y Shapley; tabla 8: las 224 combinaciones en tres especificaciones y el criterio de supervivencia) | 12 |
| `figuras/*.png`, `figuras/gris/*.png` | figuras 1 a 5, en color y en grises (las que enlaza el manuscrito) | 12 |

## 3. Códigos de causa

604 Alzheimer, Parkinson y otras del sistema nervioso; 201 tumor maligno del estómago; 613 malformaciones congénitas; 204 tumor maligno de tráquea, bronquios y pulmón; 406 otras afecciones del periodo perinatal; 108 infecciones respiratorias agudas; 210 tumor maligno de la próstata; 212 leucemia; 302 enfermedades hipertensivas; 505 ahogamiento y sumersión accidentales; 404 dificultad respiratoria del recién nacido; 103 enfermedades transmitidas por vectores y rabia; 512 agresiones (homicidios). Agrupación del DANE en la variable `causa_667` del archivo de defunciones (columna `c667` de la base extendida).
