"""Insumos propios del artículo A3 y agregados publicados (data/intermedios/): huellas, filas, columnas y totales
contra cifras publicadas por la misma fuente. Las pruebas de los microdatos del DANE están en pruebas/insumos_eevv/."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

import polars as pl
import pytest

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
from agrupamiento.constantes import CLAVES, COVS  # noqa: E402

RAW = RAIZ / "data" / "raw"
EXT = RAW / "externos"
INT = RAIZ / "data" / "intermedios"
MAN = list(csv.DictReader(open(RAW / "MANIFIESTO.csv", encoding="utf-8")))


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.mark.parametrize("f", [r for r in MAN if r["grupo"] in ("a3", "a3-config") and r["redistribuye"] == "sí"], ids=lambda r: r["archivo"])
def test_huellas_versionados(f):
    p = (RAW / f["archivo"]).resolve()
    assert p.exists(), f"falta {f['archivo']}"
    assert p.stat().st_size == int(f["bytes"])
    assert _sha(p) == f["sha256"]


def test_mgn_cnpv_personas_censadas():
    """1.122 municipios y áreas no municipalizadas; suma de personas = 44.164.417, la población censada del CNPV 2018
    (DANE, boletín técnico de la segunda entrega)."""
    g = json.loads((EXT / "63" / "mgn2018_mpio_cnpv.geojson").read_text())
    assert len(g["features"]) == 1122
    assert sum(f["properties"]["STP27_PERS"] for f in g["features"]) == 44_164_417
    fino = json.loads((EXT / "63" / "mgn2018_mpio_fino.geojson").read_text())
    assert {f["properties"]["MPIO_CDPMP"] for f in fino["features"]} == {f["properties"]["MPIO_CDPMP"] for f in g["features"]}


def test_etnia_totales_publicados():
    """Población indígena 1.905.617 y negra, mulata, afrodescendiente o afrocolombiana 2.950.072 (DANE, presentaciones
    de grupos étnicos del CNPV 2018, 2019)."""
    e = json.loads((EXT / "63" / "etnia2018_mpio.json").read_text())["features"]
    assert len(e) == 1122
    assert sum(x["attributes"]["TOTAL_UNIDADES_PA1_GRP_ETNIC_1"] or 0 for x in e) == 1_905_617
    assert sum(x["attributes"]["TOTAL_UNIDADES_PA1_GRP_ETNIC_5"] or 0 for x in e) == 2_950_072


def test_nbi_ipm():
    for f, campo in (("nbi2018_mpio.json", "NBIC_Total_Prop_de_Personas_en_NBIPorc"), ("ipm2018_mpio.json", "IPM")):
        x = json.loads((EXT / "63" / f).read_text())["features"]
        assert len(x) == 1122 and len({r["attributes"]["MPIO_CCDGO"] for r in x}) == 1122
        v = [r["attributes"][campo] for r in x]
        assert all(0 <= a <= 100 for a in v)


def test_reps_2022():
    r = pl.read_csv(EXT / "63" / "reps_capacidad_2022.csv", infer_schema=False)
    assert r.height == 41_427                      # número de filas que informa la API de datos.gov.co (count(*))
    assert r["fecha_corte"].str.contains("Nov  5 2022").all()


def test_divipola_ideam_pdet():
    assert pl.read_csv(EXT / "52" / "divipola_gdxc-w37w.csv", infer_schema=False).height == 1122
    assert pl.read_csv(EXT / "52" / "ideam_estaciones.csv", infer_schema=False).height == 9914
    cl = pl.read_csv(RAIZ / "config" / "clasificacion_municipal.csv", infer_schema=False)
    assert (cl["pdet"] == "si").sum() == 170      # Decreto Ley 893 de 2017


def test_seleccion_causas():
    s = pl.read_csv(RAIZ / "config" / "seleccion_causas_departamental.csv")
    assert s.height == 61 and (s["p_9999"] < 0.05).sum() == 27 and s["seleccionada"].sum() == 12
    assert {str(c) for c in s.filter(pl.col("seleccionada") == 1)["codigo_667"]} == {c for c in map(str, CLAVES) if c not in ("todas", "512")}


def test_intermedios_y_defunciones_publicadas():
    """Las defunciones del periodo cuadran con los totales anuales que publica el DANE, y lo que sale de la muestra está
    declarado (sin municipio, sin edad, San Andrés y Providencia)."""
    meta = json.loads((INT / "datos_meta.json").read_text())
    tot = {r["anio"]: int(r["total_publicado"]) for r in csv.DictReader(open(RAIZ / "pruebas" / "insumos_eevv" / "totales_publicados_dane.csv", encoding="utf-8"))
           if r["tipo"] == "defunciones_no_fetales"}
    for a, v in meta["perdidas"].items():
        assert v["total"] == tot[a]
    oe = pl.read_parquet(INT / "observados_esperados.parquet")
    analiz = sum(v["total"] - v["sin_mun"] - v["sin_edad"] - v["san_andres"] for v in meta["perdidas"].values())
    assert oe.filter(pl.col("causa") == "todas")["O"].sum() == analiz == 2_231_854
    for c in oe["causa"].unique():
        s = oe.filter(pl.col("causa") == c)
        assert abs(s["O"].sum() - s["E"].sum()) < 1e-6 * s["O"].sum()      # estandarización indirecta: Σ O = Σ E
    an = pl.read_parquet(INT / "observados_esperados_anual.parquet")
    assert an.group_by("causa").agg(pl.col("O").sum()).sort("causa")["O"].to_list() == oe.group_by("causa").agg(pl.col("O").sum()).sort("causa")["O"].to_list()


def test_covariables_completas():
    cov = pl.read_parquet(INT / "covariables.parquet")
    assert cov.height == 1120 and cov["mun"].n_unique() == 1120
    for v in COVS + ["basura_sin_hta", "lat", "lon", "dp", "macro"]:
        assert cov[v].null_count() == 0, v
    adj = pl.read_csv(INT / "pesos_reina.csv")
    assert adj.height == 6396 and set(adj["mun"].unique()) == set(cov["mun"])


def _cuadro5_2018():
    import pandas as pd
    p = EXT / "63" / "CUADRO5-NOFETALES-2018-def.xls"
    if not p.exists():
        pytest.skip("cuadro 5 de 2018 no descargado (scripts/00_descargar_insumos.py)")
    xl = pd.ExcelFile(p, engine="calamine")
    filas = []
    for hoja in xl.sheet_names:
        if not re.fullmatch(r"\d\d", hoja) or hoja == "00":
            continue
        x = pd.read_excel(xl, sheet_name=hoja, header=None)
        mun = None
        for i in range(len(x)):
            c0, c1, tot = x.iat[i, 0], x.iat[i, 1], x.iat[i, 2]
            if isinstance(c0, str) and re.match(r"^\d{5} ", c0):
                mun = int(c0[:5])
            elif hoja == "11" and c0 == "TOTAL":
                mun = 11001
            if mun is None or not isinstance(c1, str):
                continue
            m = re.match(r"^(\d{3}) ", c1)
            causa = "todas" if c1 == "TOTAL" else (m.group(1) if m else None)
            if causa:
                filas.append((mun, causa, int(float(tot))))
    return pl.DataFrame(filas, schema=["mun", "causa", "O_dane"], orient="row").unique(["mun", "causa"], keep="first")


def test_observados_2018_contra_cuadro5():
    """Los observados de 2018 por municipio de residencia y causa 6/67 son los del cuadro 5 del DANE, que los publica
    además por edad y sexo. Diferencias admitidas: solo las defunciones sin edad (excluidas aquí) y la fusión de Nuevo
    Belén de Bajirá (inexistente en 2018)."""
    c = _cuadro5_2018()
    an = pl.read_parquet(INT / "observados_esperados_anual.parquet").filter(pl.col("anio") == 2018).select("mun", "causa", "O")
    j = an.join(c, on=["mun", "causa"], how="left").with_columns(pl.col("O_dane").fill_null(0))
    dif = j.filter(pl.col("O") != pl.col("O_dane"))
    assert dif.height <= 3, dif.head(10)
    assert int((j["O_dane"] - j["O"]).abs().sum()) <= 3


def test_sin_registros_individuales():
    """Los intermedios son agregados municipio × causa (o × año); ninguno tiene una fila por defunción."""
    oe = pl.read_parquet(INT / "observados_esperados.parquet")
    assert oe.height == 1120 * 14 and oe.select("mun", "causa").is_unique().all()
    an = pl.read_parquet(INT / "observados_esperados_anual.parquet")
    assert an.select("anio", "mun", "causa").is_unique().all()
