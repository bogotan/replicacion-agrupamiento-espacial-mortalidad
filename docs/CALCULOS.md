# Trazabilidad de cada cifra

Mapa cifra → script y función → clave de salida → prueba. Todas las cifras del texto las comprueba `pruebas/verificar_cifras.py` (cada comprobación tiene un nombre, que se cita en la última columna); las del resumen las recalcula además, con código independiente, `pruebas/test_calculos.py`. Funciones en `src/agrupamiento/`; salidas en `output/` y `data/intermedios/`.

## Contenido

1. Ecuaciones
2. Resumen
3. Introducción
4. Datos (sección 3 y tabla 1)
5. Resultados: tablas y figuras principales
6. Conclusiones

## 1. Ecuaciones

| Ecuación | Qué es | Función | Prueba del estimador |
|---|---|---|---|
| (1) | estandarización indirecta: E = Σ λ_a P_a | `datos.py:esperados()` | `test_insumos.py:test_intermedios_y_defunciones_publicadas` (Σ O = Σ E por causa) |
| (2) | Bayes empírico global | `estimadores.py:bayes_empirico()` | `test_estimadores.py:test_bayes_empirico` |
| (3) | I de Moran | `estimadores.py:moran()`, `moran_columnas()` | `test_moran_contra_esda` |
| valor esperado bajo H0 | tr(MW)/(n−k) | `estimadores.py:esperanza_moran_h0()` | `test_esperanza_moran_h0_solo_constante` |
| (4) | valor de Shapley | `estimadores.py:shapley()` | `test_shapley_suma_y_simetria` |
| (5) | exceso sobre el placebo calibrado | `modelos.py:etapa_revision()` | `test_calculos.py` (recálculo con otra calibración y otra semilla) |
| (6), (7) | SDM y SDEM por máxima verosimilitud | `modelos.py:ajustar_causa()` (spreg `ML_Lag`, `ML_Error`) | paquete spreg 1.9.1 |
| efectos de LeSage | directo, indirecto y total | `estimadores.py:impactos_sdm()`, `impactos_sdem()` | `test_efectos_lesage_formula_espectral_contra_inversa` |
| (8) | panel con dos efectos fijos y errores agrupados | `estimadores.py:feols_dos_vias()` | `test_feols_dos_vias_contra_dummies` |
| SLX ponderado | MCP con errores agrupados (CR1) y de Conley | `estimadores.py:slx_inferencia()`, `cluster_se()`, `conley_V()` | `test_slx_errores_agrupados_contra_statsmodels` |
| BH | Benjamini-Hochberg | `estimadores.py:bh()` | `test_bh_contra_statsmodels` |
| (9) | Poisson-BYM2 escalado por Laplace | `estimadores.py:bym2_laplace()`, `escala_icar()` | `test_bym2_sin_efecto_espacial_coincide_con_poisson`, `test_bym2_recupera_parametros_simulados`, `test_escala_icar_media_geometrica_uno` |

## 2. Resumen

| Cifra en el texto | Qué es | Script:función | Clave de salida | Prueba |
|---|---|---|---|---|
| 2,23 millones de defunciones | defunciones analizadas 2018–2025 | `04_intermedios.py` → `datos.py:esperados()` | `data/intermedios/observados_esperados.parquet` (Σ O, causa `todas`) | `verificar: millones`; `test_calculos: test_tamano_muestra` |
| 1.120 municipios | municipios con vecinos terrestres | `03_pesos.py`, `04_intermedios.py` | `output/resultados.json: meta.n` | `verificar: municipios`; `test_tamano_muestra` |
| 10 % y 22 % | caída promedio de la I con 16 covariables placebo (mínimo y máximo entre las 12 series con I > 0,1) | `08_revision.py` → `modelos.py:etapa_revision()` | `revision.json: causas.<c>.caida_placebo_media` | `verificar: placebo en el resumen`; `test_placebo_entre_10_y_22` |
| entre 13 y 33 puntos más | exceso de la caída real sobre la placebo (12,5 y 32,8 redondeados) | ídem | `caida_exceso` | `verificar: exceso redondeado`; `test_exceso_13_a_33_en_11_de_12` |
| 11 de 12 series | series con p placebo ≤ 0,03 | ídem | `p_placebo` | `verificar: series con exceso`; `test_exceso_13_a_33_en_11_de_12` |
| pero no en las enfermedades hipertensivas | la excluida es 302 (exceso 5,9, p = 0,14) | ídem | `causas.302` | `verificar: la excluida es hipertensivas` |
| seis series | bloque Q con intervalo placebo por encima de cero | ídem | `bloques.Q.shapley_exceso_ic` | `verificar: Q en el resumen`; `test_calidad_del_registro_supera_al_placebo_en_seis` |
| 14 de 224 asociaciones | efecto total con q < 0,05 en SLX ponderado y en BYM2, mismo signo que el SDEM | `12_figuras_tablas.py` → `figuras.py:escribir_tablas_revision()` | `output/tablas/tabla8_efectos_coprincipales.csv: sobrevive` | `verificar: sobreviven en el resumen`; `test_sobreviven_14_de_224` |

## 3. Introducción

| Cifra | Qué es | Script:función | Clave | Prueba |
|---|---|---|---|---|
| 0,49; 1,12; 79 %; 4,28 | razón O/E por macrorregión (estómago en el Caribe y la región andina, homicidios en la Pacífica, vectores en la Orinoquía) | `figuras.py:smr_eb()` sobre `observados_esperados.parquet` | (calculada) | `verificar: estómago Caribe`, `estómago Andina`, `homicidios Pacífica`, `vectores Orinoquía` |
| tres de las 12 dejan de agruparse | p ≥ 0,05 en la I departamental recalculada | `11_departamental.py` → `modelos.py:etapa_departamental()` | `departamental.json` | `verificar: tres no se agrupan entre departamentos` |
| I de 0,12 a 0,62; 21 % y 55 %; 10 % y 22 %; 13 a 33; 14 y 13 puntos; 14 de 224 | ver resumen y sección 5 | | | `verificar: rango I 12 series`, `caída real en la introducción`, `placebo en la introducción`, `exceso redondeado intro`, `Q intro`, `sobreviven en la introducción` |

## 4. Datos (sección 3 y tabla 1)

| Cifra | Qué es | Script:función | Clave | Prueba |
|---|---|---|---|---|
| 2.260.923; 26.407 (1,2 %); 2.659 (0,12 %); 3; 2.231.854 | defunciones del periodo y exclusiones | `01_conteos.py` → `datos.py:perdidas()` | `data/intermedios/datos_meta.json` | `verificar: defunciones`, `sin municipio`, `San Andrés`, `sin edad`, `cuadre`; `test_intermedios_y_defunciones_publicadas` |
| 556 | municipios sin muertes por vectores | `06_modelos.py` | `resultados.json: causas.103.ceros` | `verificar: ceros vectores` |
| −0,52; −0,36 | correlación de la contracción con ruralidad y altitud | (calculada en la prueba) | | `verificar: C ~ ruralidad estómago`, `C ~ altitud RN` |
| 0,8 km; 93; 170; 32 valores; 0,94; 0,17; 0,09; 0,85 | estaciones, UCI, PDET, completitud, rezagos e I de las covariables | `datos.py:covariables()`, `etapa_revision()` | `covariables.parquet`, `revision.json: I_covariables` | `verificar: km estación`, `UCI`, `PDET`, `completitud valores`, `corr rezago …`, `I cov …` |
| 29,3 %; 89,3 % | completitud de Vichada y Bogotá | `02_calidad_registro.py` → `datos.py:completitud_departamental()` | `covariables.parquet: completitud` | `verificar: completitud Vichada`, `completitud Bogotá` |
| 54,8; 17,0; 32,9; 11,3 | cifras del estudio 55 (citadas, no recalculadas) | — | `config/cifras_citadas.csv` | `verificar: cifras del estudio 55` |
| tabla 1 (todas las celdas) | media, DE, P10, mediana, P90 | `figuras.py:escribir_tablas()` | `output/tablas/tabla1_covariables.csv` | `verificar: tabla 1 <variable>` |
| 5,7 vecinos; 1; 19; 11 | matriz de pesos | `pesos.py:todas()` | `resultados.json: meta.W` | `verificar: vecinos`, `forzados` |

## 5. Resultados: tablas y figuras principales

| Elemento | Qué es | Script:función | Clave | Prueba |
|---|---|---|---|---|
| Tabla 2 (14 filas) | muertes, I inicial, caída real y placebo, exceso, p, I con efectos fijos, LM robustas | `06_modelos.py`, `08_revision.py`; tabla en `figuras.py:escribir_tablas_revision()` | `revision.json: causas`, `tabla2_moran_modelos.csv` | `verificar: tabla 2 <causa>` |
| −0,007; −0,025 | E[I] bajo H0 | `etapa_revision()` | `revision.json: EI_H0` | `verificar: E[I] …` |
| 15 de 168 pasos | Shapley no monótono | ídem | `causas.<c>.no_monotono` | `verificar: pasos no monótonos` |
| 63 %, 42 %, 37 %; 59 %, 35 %, 38 % | numerador, denominador e I | ídem | `caida_num`, `caida_den`, `caida_real` | `verificar: numerador …` |
| 14,2 (11,0 a 16,1) y demás excesos de Shapley | figura 3 | ídem | `bloques.<b>.shapley_exceso`, `_ic` | `verificar: Q total`, `Q resto`, `S por causa`, `G estómago`, `S total` |
| 0,59, 0,51, 0,29 | Assunção-Reis | `ajustar_causa()` | `resultados.json: moran_ar` | `verificar: Assunção-Reis` |
| −0,09; −0,11; 0,08; 0,37; 0,31; 0,27 | I departamental recalculada y del trabajo previo | `etapa_departamental()`; `config/seleccion_causas_departamental.csv` | `departamental.json` | `verificar: I departamental recalculada`, `RN departamental`, `I departamental del trabajo previo` |
| 26 % y 50 %; 0,33 y 0,54; 22 % a 56 %; seis ajustes en la frontera | BYM2 (figura 4) | `09_bym2.py` → `etapa_bym2()` | `bym2.json` | `verificar: BYM2 …`, `frontera`, `seis ajustes en la frontera` |
| 16,7 a 151,3; 18,2 frente a 2,5; 0,44 (0,04); 0,33; 0,28; −0,03 a 0,05; 109,4 | forma de la dependencia | `ajustar_causa()`, `07_efectos_fijos.py` | `resultados.json: ols`, `fe_espacial` | `verificar: LM-R …`, `ρ con EF …`, `rango con EF` |
| Tabla 3 (panel) | δ con dos especificaciones de efectos fijos | `10_panel.py` → `etapa_panel()` | `panel.json` | `verificar: tabla 3 <causa>`, `panel …` |
| 52, 28, 14; 60 | combinaciones con q < 0,05 | `escribir_tablas_revision()` | `tabla8_efectos_coprincipales.csv` | `verificar: q < 0,05`, `q Conley` |
| Tabla 4 (13 filas) y figura 5 | efecto total en SLX ponderado y BYM2 | `figuras.py:efectos_seleccionados()` | `revision.json: slx`, `bym2.json` | `verificar: tabla 4 <causa> <variable>` |
| 24,6 % y 28,9 %; 18,9 % → 7,6 %; 10,8 % → −2,0 %; 86,3 %; −0,2 %; 43,5 % | efectos destacados | ídem | ídem | `verificar: altitud en el texto`, `códigos basura sin I10–I15`, `PDET total`, `PDET propio`, `educación y homicidios` |
| p = 0,19 | Chow de la altitud en estómago | `ajustar_causa()` | `resultados.json: regimenes` | `verificar: Chow altitud` |
| Tabla A2 | I inicial y residual por matriz | `ajustar_causa()` | `resultados.json: otras_W`; `tabla4_robustez_W.csv` | `verificar: tabla A2 <causa>` |
| 34 % a 41 %; 0,30 y 0,20; 1,33; 0,02 | robustez | `ajustar_causa()` | `otras_W`, `sin_pandemia`, `conley`, `gmm_sdem` | `verificar: homicidios entre matrices`, `sin pandemia`, `Conley de MCO`, `GMM` |
| Figuras 1 y 2 | mapas de la SMR suavizada; I inicial, residual, placebo y con efectos fijos | `figuras.py:figuras_y_tablas()` | `output/figuras/` | (gráficas de las cifras anteriores) |

## 6. Conclusiones

Las conclusiones repiten cifras ya trazadas: 11 de 12 series (`verificar: series con exceso`), la calidad del registro frente al bloque sociodemográfico en la mortalidad total (14,2 frente a 8,7; `verificar: Q total`, `S total`), la asociación circular de los códigos basura (`verificar: códigos basura sin I10–I15`) y la altitud en el cáncer de estómago (`verificar: altitud en el texto`).
