"""Arranque común de los scripts del artículo: fija la semilla de hash y un solo hilo de cálculo, y pone src/ en el
camino de importación. Sin rutas absolutas.

Determinismo (docs/REPLICACION.md, sección 5). El número de hilos de BLAS, numba y polars cambia el orden de las sumas y
con él el último dígito de los reales; con un hilo, dos corridas dan archivos idénticos byte a byte. La semilla de hash
de Python (PYTHONHASHSEED) cambia el orden en que se recorren los conjuntos; el código del artículo no depende de ella
(se comprobó corriendo el nivel A con semillas distintas), pero se fija en 0 por precaución. Si el proceso no arrancó
con estas variables, el script se vuelve a lanzar con ellas. Para la prueba de semillas distintas,
A3_SEMILLA_HASH_LIBRE=1 respeta la PYTHONHASHSEED que traiga el entorno."""
from __future__ import annotations

import os
import sys
from pathlib import Path

_HILOS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
          "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS", "POLARS_MAX_THREADS")
_libre = os.environ.get("A3_SEMILLA_HASH_LIBRE") == "1" and os.environ.get("PYTHONHASHSEED") not in (None, "", "random")
_falta = [v for v in _HILOS if os.environ.get(v) != "1"]
if _falta or (os.environ.get("PYTHONHASHSEED") != "0" and not _libre):
    for v in _HILOS:
        os.environ[v] = "1"
    if not _libre:
        os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable, *sys.argv])

_SRC = next(p / "src" for p in Path(__file__).resolve().parents if (p / "src" / "agrupamiento").is_dir())
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
