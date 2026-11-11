#!/usr/bin/env python3
"""Nivel B. Del paquete común (data/processed/ o $PROCESADO) a los conteos del artículo, en data/interim/63/:
conteos (año, municipio, grupo de edad, causa), población municipal por año y grupo de edad, proporción rural y
defunciones que quedan fuera de la muestra."""
import json

import _entorno  # noqa: F401
from agrupamiento import datos, rutas


def main():
    rutas.INTERIM.mkdir(parents=True, exist_ok=True)
    datos.leer_defunciones().write_parquet(rutas.INTERIM / "conteos.parquet")
    pob, rur = datos.leer_poblacion()
    pob.write_parquet(rutas.INTERIM / "poblacion.parquet"); rur.write_parquet(rutas.INTERIM / "rural.parquet")
    (rutas.INTERIM / "perdidas.json").write_text(json.dumps(datos.perdidas(), indent=1))
    datos.departamental_2022().write_parquet(rutas.INTERIM / "departamental_2022.parquet")
    print("  conteos, población, pérdidas y departamental 2022 listos")


if __name__ == "__main__":
    main()
