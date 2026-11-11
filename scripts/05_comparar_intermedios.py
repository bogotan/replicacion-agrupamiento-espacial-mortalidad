#!/usr/bin/env python3
"""Nivel B. Compara los intermedios reconstruidos (--nuevo) con los versionados (--referencia): mismas columnas, filas
y valores; los números en coma flotante se aceptan con diferencia relativa menor de 1e-9."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl


def igual(a: pl.DataFrame, b: pl.DataFrame) -> str:
    if a.columns != b.columns:
        return f"columnas distintas: {a.columns} / {b.columns}"
    if a.height != b.height:
        return f"filas distintas: {a.height} / {b.height}"
    for c in a.columns:
        x, y = a[c], b[c]
        if x.dtype.is_float():
            xa, ya = x.to_numpy(), y.to_numpy()
            ok = np.isclose(xa, ya, rtol=1e-9, atol=1e-12, equal_nan=True)
            if not ok.all():
                return f"{c}: {int((~ok).sum())} valores distintos (máx. {np.nanmax(np.abs(xa - ya)):.3g})"
        elif not x.equals(y):
            return f"{c}: valores distintos"
    return "idéntico"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nuevo", required=True); ap.add_argument("--referencia", required=True)
    a = ap.parse_args(); N, R = Path(a.nuevo), Path(a.referencia)
    fallas = 0
    for f in sorted(R.iterdir()):
        n = N / f.name
        if not n.exists():
            print(f"  FALTA    {f.name}"); fallas += 1; continue
        if f.suffix == ".parquet":
            r = igual(pl.read_parquet(n), pl.read_parquet(f))
        elif f.suffix == ".csv":
            r = igual(pl.read_csv(n), pl.read_csv(f))
        else:
            r = "idéntico" if json.loads(n.read_text()) == json.loads(f.read_text()) else "contenido distinto"
        print(f"  {'ok   ' if r == 'idéntico' else 'FALLA'}    {f.name}: {r}")
        fallas += r != "idéntico"
    print("Intermedios idénticos." if not fallas else f"{fallas} archivos distintos.")
    sys.exit(1 if fallas else 0)


if __name__ == "__main__":
    main()
