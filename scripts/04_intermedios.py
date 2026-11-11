#!/usr/bin/env python3
"""Nivel B. Escribe los agregados publicables que alimentan el nivel A (por defecto en data/intermedios/; con --salida,
en otra carpeta para compararlos con los versionados):

  observados_esperados.parquet                muertes observadas y esperadas por municipio y causa, 2018–2025
  observados_esperados_sin_pandemia.parquet   ídem sin 2020 y 2021
  observados_esperados_anual.parquet          ídem por año (panel)
  covariables.parquet                         covariables municipales y coordenadas
  departamental_2022.parquet                  observados y esperados por departamento, 2022
  pesos_reina.csv, pesos_departamentos.csv    listas de adyacencia
  datos_meta.json                             defunciones fuera de la muestra y número de municipios
"""
import argparse
import json
import shutil
from pathlib import Path

import polars as pl

import _entorno  # noqa: F401
from agrupamiento import datos, rutas
from agrupamiento.constantes import ANIOS, EXCLUIR, SIN_PANDEMIA


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--salida", default=str(rutas.INTER))
    sal = Path(ap.parse_args().salida); sal.mkdir(parents=True, exist_ok=True)
    I = rutas.INTERIM
    conteos = pl.read_parquet(I / "conteos.parquet"); pob = pl.read_parquet(I / "poblacion.parquet")
    muns = sorted(set(pob["mun"].unique().to_list()) - EXCLUIR)
    datos.esperados(conteos, pob, muns, ANIOS).write_parquet(sal / "observados_esperados.parquet")
    datos.esperados(conteos, pob, muns, SIN_PANDEMIA).write_parquet(sal / "observados_esperados_sin_pandemia.parquet")
    pl.concat([datos.esperados(conteos, pob, muns, [a]).with_columns(pl.lit(a).alias("anio")) for a in ANIOS]) \
        .select("anio", "mun", "causa", "O", "E").write_parquet(sal / "observados_esperados_anual.parquet")
    cov = datos.covariables(muns, pl.read_parquet(I / "rural.parquet"), pob, pl.read_parquet(I / "calidad.parquet"),
                            pl.read_parquet(I / "completitud.parquet"), pl.read_parquet(I / "basura_sin_hta.parquet"))
    cov.write_parquet(sal / "covariables.parquet")
    shutil.copy2(I / "departamental_2022.parquet", sal / "departamental_2022.parquet")
    for f in ("pesos_reina.csv", "pesos_departamentos.csv"):
        shutil.copy2(I / f, sal / f)
    (sal / "datos_meta.json").write_text(json.dumps({"perdidas": json.loads((I / "perdidas.json").read_text()), "n_mun": len(muns)}, indent=1))
    print(f"  intermedios en {sal}")


if __name__ == "__main__":
    main()
