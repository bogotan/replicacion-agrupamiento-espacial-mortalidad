#!/usr/bin/env python3
"""Nivel B. Calidad del registro: indicadores municipales 2018–2025 (definiciones del estudio 55), completitud
departamental 2015–2017 frente a las estimaciones de mortalidad del DANE y códigos basura sin I10–I15."""
import polars as pl

import _entorno  # noqa: F401
from agrupamiento import datos, rutas


def main():
    rutas.INTERIM.mkdir(parents=True, exist_ok=True)
    datos.calidad_municipal().write_parquet(rutas.INTERIM / "calidad.parquet")
    datos.completitud_departamental().write_parquet(rutas.INTERIM / "completitud.parquet")
    muns = pl.read_parquet(rutas.INTERIM / "rural.parquet")["mun"].sort().to_list()
    datos.basura_sin_hipertension(muns).write_parquet(rutas.INTERIM / "basura_sin_hta.parquet")
    print("  calidad del registro lista")


if __name__ == "__main__":
    main()
