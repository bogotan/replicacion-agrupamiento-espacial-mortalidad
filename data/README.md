# data/

| Carpeta | Qué hay | ¿Versionado? | Cómo se produce |
|---|---|---|---|
| `raw/MANIFIESTO.csv` | una fila por insumo: fuente, URL, fecha, tamaño, SHA-256, licencia, si se redistribuye y quién lo usa | sí | a mano; `scripts/00_descargar_insumos.py` lo verifica |
| `raw/externos/63/` | CNPV 2018 y MGN 2018 (DANE), NBI, IPM, etnia, REPS 2022; el cuadro 5 de defunciones 2018 se descarga | sí, salvo el cuadro 5 | descarga documentada en `raw/externos/63/FUENTE.md` |
| `raw/externos/52/` | DIVIPOLA con coordenadas y catálogo de estaciones del IDEAM (del estudio 52) | sí | `raw/externos/52/FUENTE.md` |
| `raw/*.zip`, `raw/antiguos/`, `raw/poblacion/` | microdatos y proyecciones del DANE | no | `scripts/00_descargar_insumos.py` (llama a `scripts/insumos_eevv/00_descargar_microdatos.py`) |
| `processed/` | bases del paquete común (una fila por evento) | no | `make bases-eevv` |
| `interim/63/` | conteos por edad, población, calidad del registro, pesos (nivel B) | no | `make intermedios` |
| `intermedios/` | agregados publicables que alimentan el nivel A | sí | `make intermedios NUEVOS=data/intermedios` (o `scripts/04_intermedios.py`) |

Variables: `docs/DICCIONARIO.md`. Procedencia y términos: `docs/INSUMOS.md`.
