# ¿Cuánto del agrupamiento espacial de la mortalidad por causa en Colombia explican los determinantes y la calidad del registro?

Repositorio de replicación del artículo A3. Autor: Daniel A. Garavito Jiménez (Sinergística, Bogotá). Octubre de 2026.

**Estado.** Artículo de Daniel A. Garavito Jiménez (Sinergística S.A.S., ORCID 0009-0009-6871-5807), enviado a *Cuadernos de Economía* (Universidad Nacional de Colombia) en noviembre de 2026. El texto del artículo y sus entregables se agregan a este repositorio cuando se publique; mientras tanto, el código, los agregados y las pruebas reproducen todas las cifras. En la preparación se usó Claude (Anthropic) como asistente de código y redacción, bajo la dirección y revisión del autor; los registros de `docs/REVISION.md` son revisiones simuladas por esa misma herramienta, no arbitraje externo.

## Contribuciones

**Daniel A. Garavito Jiménez** (autor):

- definió la pregunta del artículo y eligió la revista;
- escogió las fuentes y aprobó el diseño: unidad de análisis, tratamiento o escenarios, especificaciones e inferencia;
- tomó las decisiones de privacidad y de licencia: no publicar celdas de 1 a 4 en cruces sensibles y licenciar por capa (código, método, datos y texto);
- aportó su conocimiento del sistema de salud y de las estadísticas públicas de Colombia, buena parte de él adquirido como asesor de Innovación y Analítica de la ADRES;
- revisó los resultados, ordenó las correcciones y responde por todo el contenido.

**Claude (Anthropic)**, herramienta de inteligencia artificial usada bajo su dirección: escribió y depuró el código, ejecutó los análisis, produjo cuadros y figuras, redactó borradores de la documentación y simuló revisiones de árbitro. No figura como autor.

## Contenido

1. La pregunta y el hallazgo
2. Cómo reproducir
3. Tiempos y memoria
4. Mapa del repositorio
5. Insumos, licencias y privacidad
6. Limitaciones conocidas

## 1. La pregunta y el hallazgo

¿Por qué se parecen los municipios vecinos en su mortalidad por causa? El artículo mide, para la mortalidad total, los homicidios y 12 causas de la lista 6/67, cuánto de la autocorrelación espacial municipal (2.231.854 defunciones de 2018 a 2025, 1.120 municipios) desaparece al condicionar en determinantes sociodemográficos, geografía, oferta de servicios y calidad del registro. Covariables placebo con la misma estructura espacial reducen la I de Moran entre 10 % y 22 %; las reales la reducen entre 13 y 33 puntos más en 11 de 12 series, pero no en las enfermedades hipertensivas. La calidad del registro aporta más que el placebo en seis series. Con estos datos no se distingue el contagio de los determinantes compartidos. De 224 asociaciones entre causa y covariable, 14 sobreviven en las dos especificaciones coprincipales (SLX ponderado con errores agrupados y Poisson con efecto BYM2).

## 2. Cómo reproducir

Requisitos: Python 3.13 con `pip install -r requirements.txt`; para los entregables, pandoc y XeLaTeX. Todo corre desde la raíz con `make`.

**Nivel A, desde los agregados versionados (`data/intermedios/`).** Reproduce todas las estimaciones, tablas, figuras y cifras del texto.

```bash
make nivel-a          # modelos (unos 22 minutos, la mayor parte el BYM2), figuras, tablas y verificación de cifras
make pruebas          # pytest (insumos versionados, estimadores, recálculo independiente del resumen) y verificación
make nivel-a-rapido   # sin reestimar: usa output/*.json versionados; figuras, tablas y verificación (menos de 1 minuto)
```

`pruebas/verificar_cifras.py` termina en «Todo coincide.» si cada cifra del manuscrito coincide con las salidas (requiere el texto del artículo, que se agrega cuando se publique; sin él, el verificador lo avisa y termina).

**Nivel B, desde los insumos originales.** Descarga y verifica los microdatos del DANE, reconstruye las bases con el paquete común y los agregados del artículo, y los compara con los versionados.

```bash
make insumos          # baja lo que falte y verifica huellas SHA-256 (unos 3 GB de zips del DANE)
make bases-eevv       # paquete común: data/processed/ (unos 10 minutos)
make pruebas-eevv     # huellas, columnas y totales del DANE
make intermedios      # etapas 01 a 04 -> data/interim/63/intermedios_b/
make comparar         # mismas filas y valores que data/intermedios/
```

O `make nivel-b`. Si ya se tienen las bases del paquete común en otra carpeta, basta con `PROCESADO=<carpeta processed> make intermedios comparar` (y los archivos de `data/raw/poblacion/dane_mor/`).

**Etapas.**

| Script | Nivel | Qué hace |
|---|---|---|
| `00_descargar_insumos.py` | B | descarga y verifica huellas (grupos «a3», «a3-config» y «eevv» del manifiesto) |
| `insumos_eevv/01…04` | B | paquete común: zips del DANE → `data/processed/` |
| `01_conteos.py` | B | conteos por municipio, edad y causa; población; pérdidas; departamental 2022 |
| `02_calidad_registro.py` | B | calidad del registro municipal, completitud departamental, códigos basura sin I10–I15 |
| `03_pesos.py` | B | contigüidad tipo reina desde el MGN 2018 y contigüidad departamental |
| `04_intermedios.py` | B | escribe `data/intermedios/` (o `--salida`) |
| `05_comparar_intermedios.py` | B | compara intermedios reconstruidos con los versionados |
| `06_modelos.py` | A | Moran, MCO por bloques, Shapley, LM, SAR, SEM, SDM, SDEM, SLX, efectos, Conley, GMM, regímenes |
| `07_efectos_fijos.py` | A | SDM y SDEM con efectos fijos de departamento |
| `08_revision.py` | A | placebo calibrado, E[I], numerador y denominador, SLX ponderado con errores agrupados y de Conley, BH |
| `09_bym2.py` | A | Poisson-BYM2 escalado con rezagos de las covariables (Laplace) |
| `10_panel.py` | A | panel anual con rezago espacio-temporal |
| `11_departamental.py` | A | I de Moran departamental de 2022 (comparación con el trabajo previo) |
| `12_figuras_tablas.py` | A | figuras (color y grises) y tablas |

Semillas: 20261024 (efectos del SDM y del SAR, permutaciones) y 20261025 (placebo). Detalle en `docs/REPLICACION.md`.

## 3. Tiempos y memoria

Corridas en copias limpias del repositorio, 30 de octubre de 2026, máquina de 2 núcleos compartida (detalle en `docs/REPLICACION.md`).

| Corrida | Duración | Memoria máxima |
|---|---|---|
| Nivel A completo (`make nivel-a`) | 1.289 s (21,5 min; el BYM2, 1.122 s) | 536 MB |
| Pruebas (`make pruebas`) | 64 s | 606 MB |
| Nivel A rápido (`make nivel-a-rapido`) | 8 s | 475 MB |
| Nivel B, etapas propias (`make intermedios comparar`) | unos 8 s | 573 MB |
| Paquete común, desde los zips | 604 s | 932 MB (`docs/REPLICACION_EEVV.md`) |

## 4. Mapa del repositorio

```
├── README.md, LICENSE (MIT), LICENSE-datos.md (CC BY-NC 4.0, CC BY 4.0 para la documentación de método y términos de las fuentes), CITATION.cff, requirements.txt, Makefile
├── config/          clasificación PDET, códigos basura, selección de causas del trabajo previo, cifras citadas
├── src/agrupamiento/   rutas, constantes, datos (nivel B), pesos, estimadores, modelos (nivel A), figuras y tablas
├── src/eevv_insumos/   paquete común de Estadísticas Vitales (sin cambios)
├── scripts/         00–12 del artículo; insumos_eevv/ del paquete común
├── data/            raw/ (MANIFIESTO.csv y externos abiertos), intermedios/ (versionados), processed/ e interim/ (no)
├── output/          resultados.json, revision.json, bym2.json, panel.json, departamental.json, tablas/, figuras/ (y gris/)
├── pruebas/         test_insumos.py, test_estimadores.py, test_calculos.py, verificar_cifras.py, insumos_eevv/
└── docs/            INSUMOS.md, CALCULOS.md, DICCIONARIO.md, REPLICACION.md, REVISION.md, FUENTES_LITERATURA.md,
                     y los documentos del paquete común (INSUMOS_EEVV.md, DICCIONARIO_EEVV.md, REPLICACION_EEVV.md)
```

## 5. Insumos, licencias y privacidad

- Microdatos del DANE: no se redistribuyen (términos del Archivo Nacional de Datos); se descargan y verifican por su huella.
- Insumos abiertos pequeños (MGN y CNPV 2018, NBI, IPM, etnia, REPS 2022, DIVIPOLA, IDEAM, PDET): versionados en `data/raw/externos/` con su ficha en `docs/INSUMOS.md`.
- `data/intermedios/` tiene solo agregados municipio × causa (y × año) y covariables municipales, sin registros individuales. Los observados por municipio y causa tienen la misma o menor granularidad que el cuadro 5 del DANE; la decisión sobre las proporciones de calidad del registro con numeradores de 1 a 4 está en `docs/INSUMOS.md`, sección 4.
- Código: MIT. Texto, tablas, figuras y datos derivados: CC BY-NC 4.0 (con CC BY-SA 4.0 para las columnas derivadas de datos.gov.co). Documentación de método: CC BY 4.0.

## 6. Limitaciones conocidas

- Proporciones de calidad del registro con numeradores de 1 a 4 en 42 municipios: se publican sin supresión porque son proporciones sin conteos y no describen a personas (`docs/INSUMOS.md`, sección 4).
- El BYM2 se estima por Laplace (sin INLA ni MCMC); seis ajustes quedan con φ en la frontera.
