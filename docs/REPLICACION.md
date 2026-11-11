# Replicación

Cómo se reproduce el artículo A3 y qué resultó al hacerlo en copias limpias del repositorio. Los insumos y sus huellas están en `docs/INSUMOS.md`; los cálculos, en `docs/CALCULOS.md`.

## Contenido

1. Resultado de las corridas en limpio
2. Máquina y versiones
3. Nivel A
4. Nivel B
5. Semillas y aleatoriedad
6. Cómo repetirlo
7. Diferencias con el código anterior

## 1. Resultado de las corridas en limpio

| Corrida | Fecha (UTC) | Resultado |
|---|---|---|
| Nivel A (`make nivel-a`) | 30 de octubre de 2026, 00:11 a 00:32 | «201 comprobaciones. Todo coincide.»; JSON y tablas idénticos byte a byte a los versionados de entonces (en `bym2.json` cambia solo el campo de segundos) |
| Pruebas (`make pruebas`) | 30 de octubre de 2026, 00:35 | 44 pruebas pasan, 1 omitida (la del cuadro 5 cuando no está descargado); «Todo coincide.». Tras la sección 5: 45 pruebas pasan con el cuadro 5 descargado, 202 comprobaciones |
| Nivel A rápido (`make nivel-a-rapido`) | 30 de octubre de 2026 | 8 s; «Todo coincide.» |
| Nivel A dos veces, con `PYTHONHASHSEED` 0 y 4242 (`make nivel-a`) | 30 de octubre de 2026, 02:32 a 03:03 | salidas no gráficas idénticas byte a byte entre las dos corridas, salvo los segundos de cómputo de `bym2.json`; figuras idénticas; «202 comprobaciones. Todo coincide.» (sección 5) |
| Nivel B dos veces, con `PYTHONHASHSEED` 0 y 4242 (`make intermedios`) | 30 de octubre de 2026, 02:32 | los 8 intermedios idénticos byte a byte entre sí y con los versionados |
| Nivel B, etapas propias (`make intermedios comparar`) | 30 de octubre de 2026, 00:36 a 00:37 | 20 huellas verificadas; cuadro 5 descargado con huella correcta; los 8 archivos de `data/intermedios/` reconstruidos «idéntico» (enteros iguales; reales con tolerancia relativa 1e-9, diferencia máxima observada 1e-11 en E); `test_insumos.py`: 28 pruebas pasan |
| Paquete común desde los zips del DANE | 29 de octubre de 2026 | 604 s, 932 MB; ver `docs/REPLICACION_EEVV.md` |

Cada corrida se hizo en una copia del repositorio sin `output/` regenerable ni `data/interim/`, con `nice -n 10`.

## 2. Máquina y versiones

2 núcleos Intel Xeon a 2,1 GHz compartidos con otros cuatro procesos, 7,8 GB de RAM, Linux, Python 3.13. Versiones exactas en `requirements.txt` (spreg 1.9.1, libpysal 4.15.1, esda 2.10.0, numpy 2.5.3, scipy 1.18.1, statsmodels, polars 2.0.0, geopandas 1.2.0, matplotlib 3.11.2, entre otras).

## 3. Nivel A

Desde `data/intermedios/` (versionado). Tiempos y memoria máxima por etapa (`tiempos_nivel-a.tsv`, escrito por `scripts/insumos_eevv/medir.py`):

| Etapa | Segundos | Memoria máxima (MB) |
|---|--:|--:|
| `06_modelos.py` | 108 | 460 |
| `07_efectos_fijos.py` | 27 | 415 |
| `08_revision.py` | 21 | 536 |
| `09_bym2.py` | 1.122 | 412 |
| `10_panel.py` | 3 | 360 |
| `11_departamental.py` | 0 | 83 |
| `12_figuras_tablas.py` | 6 | 475 |
| Total | 1.289 (21,5 min) | 536 |
| `make pruebas` | 64 | 606 |

El BYM2 ocupa el 87 % del tiempo (14 series, modelo nulo y completo, varios arranques de Nelder-Mead por ajuste).

## 4. Nivel B

Desde los insumos originales. Con las bases del paquete común ya construidas (`PROCESADO=<carpeta>/data/processed`):

| Etapa | Segundos | Memoria máxima (MB) |
|---|--:|--:|
| `00_descargar_insumos.py --sin-eevv` (verificación y cuadro 5) | 2 | 31 |
| `01_conteos.py` | 1 | 573 |
| `02_calidad_registro.py` | 1 | 194 |
| `03_pesos.py` | 3 | 479 |
| `04_intermedios.py` | 1 | 163 |

`05_comparar_intermedios.py` compara fila por fila con `data/intermedios/`: covariables y observados idénticos; esperados con diferencia relativa máxima de 1e-11; `datos_meta.json` igual.

## 5. Semillas, hilos y determinismo

- 20261024: extracciones de los efectos del SDM y del SAR (2.000; 1.000 en las matrices alternativas) y permutaciones de la I de Moran (999), fijadas con `np.random.seed` antes de cada llamada a `esda`.
- 20261025: placebo (calibración de ρ y 200 réplicas).
- `pruebas/test_calculos.py` recalcula el placebo con semilla 99, calibración por bisección y 100 réplicas, sin usar el paquete del artículo, y exige las mismas conclusiones con 1,5 puntos de tolerancia.

Los valores p por permutación tienen resolución de 0,001.

**Protecciones.** En A1, A2 y A4 aparecieron diferencias entre corridas por el orden de recorrido de conjuntos (semilla de hash), por los hilos de BLAS y numba y por las sumas paralelas de `group_by` en polars. Aquí:

- el Makefile exporta `PYTHONHASHSEED=0` y un hilo en `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `NUMEXPR_NUM_THREADS`, `NUMBA_NUM_THREADS` y `POLARS_MAX_THREADS`; `scripts/_entorno.py` relanza el script con esas variables si no las trae (`A3_SEMILLA_HASH_LIBRE=1` respeta otra semilla de hash, solo para la prueba de abajo);
- las sumas de reales por grupo (esperados municipales y departamentales, esperados de la completitud) se hacen con numpy después de ordenar (`datos.suma_ordenada`); las demás sumas son de enteros (conteos y población, que es entera) y son exactas en cualquier orden;
- las uniones de polars declaran `maintain_order="left"` y el panel se ordena por municipio y año;
- un ajuste del panel sin variación dentro de los efectos fijos (la leucemia) se declara no estimable en vez de informar el ruido de redondeo, que cambiaba con el último dígito de los esperados.

**Prueba.** En dos copias limpias se corrieron `make intermedios` y `make nivel-a` completos en paralelo, una con `PYTHONHASHSEED=0` y otra con 4242. Resultado: los 8 intermedios y todas las salidas no gráficas (`resultados.json`, `revision.json`, `panel.json`, `departamental.json`, las 9 tablas) son idénticos byte a byte; `bym2.json` difiere solo en el campo `seg` (segundos de cómputo); las 10 figuras son idénticas byte a byte. Tiempos por etapa similares a los de la sección 3 (06: 108 s; 09: 1.164 s), memoria máxima 536 MB por corrida. La primera pasada se detuvo en el panel por el caso de la leucemia; tras la corrección, las etapas 10 a 12 y la verificación se volvieron a correr en las dos copias (las etapas 06 a 09 no usan el código corregido).

**Tolerancia frente a la versión anterior.** No hay tolerancia entre corridas con el mismo código. Al pasar del commit inicial a este, los esperados cambiaron 1 o 2 unidades del último dígito (orden de suma) y BLAS pasó de dos hilos a uno. Eso movió las salidas así, sin cambiar ninguna cifra del texto ni un píxel de las figuras:

| Salida | Diferencia absoluta máxima | Origen |
|---|---|---|
| `revision.json`, `panel.json`, `departamental.json` | 2e-14 | redondeo |
| `resultados.json` | 3e-7 (valores p de los efectos del SDEM con matrices alternativas) | redondeo amplificado por las simulaciones de efectos |
| `bym2.json`, β | 2e-8 | ídem, en los ajustes con φ en la frontera |
| `bym2.json`, θ (logit de φ) | 0,42 en el cáncer de pulmón (completo); 0,10 en las enfermedades por vectores | verosimilitud marginal plana con φ ≥ 0,998 (frontera): Nelder-Mead se detiene en otro punto de la meseta; φ, σ² y β iguales a cuatro decimales |

## 6. Cómo repetirlo

```bash
pip install -r requirements.txt
make nivel-a && make pruebas            # desde los agregados
make nivel-b                            # desde los insumos (descarga unos 3 GB de zips del DANE)
PROCESADO=/ruta/data/processed make intermedios comparar   # si ya se tienen las bases comunes
```

`make limpiar` borra lo regenerable. `python scripts/09_bym2.py --rehacer` reestima el BYM2 aunque exista `output/bym2.json`.

## 7. Diferencias con el código anterior

El código venía de un solo script (`63_econometria_espacial.py`, en el repositorio de Estadísticas Vitales). Al reestructurarlo:

- se partió en las etapas 01 a 12, con las funciones en `src/agrupamiento/` y pruebas contra casos conocidos;
- los insumos de los estudios 52 y 55 se producen dentro del repositorio (`02_calidad_registro.py`) o se versionan con su origen (`config/`);
- se quitó el modelo de Leroux sin escalar de la versión 1 (lo reemplaza el BYM2) y su tabla;
- las permutaciones de `esda` quedaron con semilla fija; el valor p de las malformaciones pasó de 0,002 a 0,001 y el manuscrito lo refleja;
- se agregó la I de Moran departamental de 2022 (`11_departamental.py`).

Los cambios de la revisión de árbitro están en `docs/REVISION.md`.
