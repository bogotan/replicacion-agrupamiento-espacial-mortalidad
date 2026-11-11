#!/usr/bin/env python3
"""Nivel B. Matriz de contigüidad tipo reina desde los polígonos del MGN 2018 y contigüidad departamental derivada."""
import polars as pl

import _entorno  # noqa: F401
from agrupamiento import pesos, rutas
from agrupamiento.constantes import EXCLUIR


def main():
    muns = sorted(set(pl.read_parquet(rutas.INTERIM / "poblacion.parquet")["mun"].unique().to_list()) - EXCLUIR)
    adj = pesos.reina_desde_mgn(muns)
    adj.write_csv(rutas.INTERIM / "pesos_reina.csv")
    pesos.departamentos_desde_municipios(adj).write_csv(rutas.INTERIM / "pesos_departamentos.csv")
    print(f"  {len(muns)} municipios, {adj.height} pares de vecinos")


if __name__ == "__main__":
    main()
