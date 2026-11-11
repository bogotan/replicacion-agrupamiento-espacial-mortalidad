"""Raíz y carpetas del repositorio. Sin rutas absolutas.

La raíz es la carpeta que contiene `src/agrupamiento/`, o la que indique la variable de entorno `RAIZ`.
`PROCESADO` (variable de entorno) permite leer las bases del paquete común desde otra carpeta en el nivel B.
"""
from __future__ import annotations

import os
from pathlib import Path


def _raiz() -> Path:
    if os.environ.get("RAIZ"):
        return Path(os.environ["RAIZ"]).resolve()
    return next(p for p in Path(__file__).resolve().parents if (p / "src" / "agrupamiento").is_dir())


RAIZ = _raiz()
CONFIG = RAIZ / "config"
RAW = RAIZ / "data" / "raw"
EXT63 = RAW / "externos" / "63"
EXT52 = RAW / "externos" / "52"
PROCESADO = Path(os.environ["PROCESADO"]).resolve() if os.environ.get("PROCESADO") else RAIZ / "data" / "processed"
INTERIM = RAIZ / "data" / "interim" / "63"
INTER = RAIZ / "data" / "intermedios"
OUT = RAIZ / "output"
TAB = OUT / "tablas"
FIG = OUT / "figuras"
