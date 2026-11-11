"""Del paquete común (data/processed) y de los insumos externos a los agregados del artículo (nivel B).

Funciones:
  leer_defunciones()          conteos (año, municipio, grupo de edad, causa) 2018–2025 por residencia
  leer_poblacion()            población municipal por año y grupo de edad; proporción rural
  esperados()                 muertes observadas y esperadas por estandarización indirecta (ecuación 1)
  perdidas()                  defunciones que quedan fuera de la muestra
  calidad_municipal()         indicadores municipales de calidad del registro (como el estudio 55)
  completitud_departamental() completitud del registro 2015–2017 frente a las estimaciones del DANE (como el script 35
                              del repositorio EEVV, del que viene el estudio 55)
  basura_sin_hipertension()   códigos basura sin las defunciones I10–I15
  departamental_2022()        observados y esperados por departamento en 2022 (comparación con el trabajo previo)
  covariables()               covariables municipales de fuentes abiertas
"""
from __future__ import annotations

import json
import math

import numpy as np
import polars as pl

from . import rutas
from .codigos import clasificador
from .constantes import ANIOS, CLAVES, EXCLUIR, FUSION, MACRO, clave

G = 19  # grupos de edad: 0, 1-4, 5-9, …, 80-84, 85+


def suma_ordenada(df: pl.DataFrame, por: str, valor: pl.Expr, nombre: str) -> pl.DataFrame:
    """Suma de reales por grupo en un orden fijo: ordena por (grupo, valor) y suma con numpy. Las sumas de group_by de
    polars se reparten entre hilos y su orden cambia entre corridas, lo que mueve el último dígito (docs/REPLICACION.md)."""
    t = df.select(pl.col(por), valor.cast(pl.Float64).alias("_v")).sort([por, "_v"])
    k, v = t[por].to_numpy(), t["_v"].to_numpy()
    ini = np.flatnonzero(np.r_[True, k[1:] != k[:-1]]) if len(k) else np.array([], dtype=int)
    return pl.DataFrame({por: k[ini], nombre: np.add.reduceat(v, ini) if len(k) else v}).with_columns(
        pl.col(por).cast(df.schema[por]))


def mun_fusion(c: pl.Expr) -> pl.Expr:
    return c.replace(FUSION)


def grupo_edad(e: pl.Expr) -> pl.Expr:
    return pl.when(e == 0).then(0).when(e < 5).then(1).otherwise(pl.min_horizontal(e // 5 + 1, pl.lit(18))).cast(pl.Int8)


def _def(anio: int, cols: list[str]) -> pl.LazyFrame:
    return pl.scan_parquet(rutas.PROCESADO / "ext" / f"def_{anio}.parquet").select(cols)


def leer_defunciones(anios=ANIOS) -> pl.DataFrame:
    """Conteos (año, municipio, grupo de edad, causa) de las 14 series; residencia; sin municipio o sin edad se excluyen."""
    cods = [c for c in CLAVES if c != "todas"]
    partes = []
    for a in anios:
        d = (_def(a, ["mun", "g", "c667"]).filter((pl.col("mun") > 0) & (pl.col("g") >= 0))
             .with_columns(mun_fusion(pl.col("mun")).alias("mun"), pl.min_horizontal(pl.col("g"), pl.lit(18)).cast(pl.Int8).alias("g"))
             .collect())
        tot = d.group_by(["mun", "g"]).agg(pl.len().alias("o")).with_columns(pl.lit("todas").alias("causa"))
        esp = (d.filter(pl.col("c667").is_in(cods)).group_by(["mun", "g", "c667"]).agg(pl.len().alias("o"))
               .with_columns(pl.col("c667").cast(pl.Utf8).alias("causa")).drop("c667"))
        partes.append(pl.concat([tot, esp.select(tot.columns)]).with_columns(pl.lit(a).alias("anio")))
    return pl.concat(partes).sort(["anio", "mun", "g", "causa"])


def leer_poblacion(anios=ANIOS) -> tuple[pl.DataFrame, pl.DataFrame]:
    p = (pl.scan_parquet(rutas.PROCESADO / "mun_pob_edad.parquet")
         .filter(pl.col("ano").is_between(anios[0], anios[-1]))
         .with_columns(mun_fusion(pl.col("mpio").cast(pl.Int32)).alias("mun"), grupo_edad(pl.col("edad")).alias("g"))
         .collect())
    tot = p.filter(pl.col("area") == "Total").group_by(["ano", "mun", "g"]).agg(pl.col("pob").sum())
    rur = (p.filter(pl.col("area") != "Total").group_by(["mun"])
           .agg((pl.col("pob").filter(pl.col("area").str.starts_with("Centros")).sum() / pl.col("pob").sum()).alias("rural")))
    return tot.rename({"ano": "anio"}), rur


def esperados(conteos: pl.DataFrame, pob: pl.DataFrame, muns: list[int], anios: list[int]) -> pl.DataFrame:
    """Ecuación (1): E_ic = Σ_a λ_ac P_ia, con λ_ac la tasa de la muestra (Σ_j O_jac / Σ_j P_ja)."""
    c = conteos.filter(pl.col("anio").is_in(anios) & pl.col("mun").is_in(muns)).group_by(["mun", "g", "causa"]).agg(pl.col("o").sum())
    p = pob.filter(pl.col("anio").is_in(anios) & pl.col("mun").is_in(muns)).group_by(["mun", "g"]).agg(pl.col("pob").sum())
    filas = []
    for k in CLAVES:
        ck = c.filter(pl.col("causa") == clave(k))
        tasa = (ck.group_by("g").agg(pl.col("o").sum()).join(p.group_by("g").agg(pl.col("pob").sum()), on="g", how="right")
                .with_columns((pl.col("o").fill_null(0) / pl.col("pob")).alias("lam")).select(["g", "lam"]))
        e = suma_ordenada(p.join(tasa, on="g"), "mun", pl.col("lam") * pl.col("pob"), "E")
        o = ck.group_by("mun").agg(pl.col("o").sum().alias("O"))
        filas.append(pl.DataFrame({"mun": muns}).join(e, on="mun", how="left", maintain_order="left").join(o, on="mun", how="left", maintain_order="left")
                     .with_columns(pl.col("O").fill_null(0), pl.col("E").fill_null(0.0), pl.lit(clave(k)).alias("causa")))
    return pl.concat(filas).select("mun", "causa", "O", "E")


def perdidas(anios=ANIOS) -> dict:
    out = {}
    for a in anios:
        d = _def(a, ["mun", "g"]).collect()
        out[str(a)] = {"total": d.height, "sin_mun": int((d["mun"] == 0).sum()), "sin_edad": int(((d["mun"] > 0) & (d["g"] < 0)).sum()),
                       "san_andres": int(d["mun"].is_in(list(EXCLUIR)).sum())}
    return out


def _clasificador_basura():
    B = pl.read_csv(rutas.CONFIG / "codigos_basura_naghavi2010.csv")
    return clasificador([(r["tipo"], r["cie10"], "K71.7" if r["tipo"] == 2 else "") for r in B.iter_rows(named=True)])


def calidad_municipal(anios=ANIOS) -> pl.DataFrame:
    """Indicadores municipales de calidad del registro, 2018–2025, con las definiciones del estudio 55:
    sin_asistencia = muertes por causas no externas sin asistencia médica / muertes por causas no externas;
    basura = muertes con código basura (Naghavi et al., 2010, tipos 1–4) / muertes con causa CIE-10;
    casa = muertes en casa o domicilio / todas las muertes. Residencia; municipio desconocido excluido."""
    fb = _clasificador_basura()
    partes = []
    for a in anios:
        d = _def(a, ["mun", "cb", "sit_defun", "asis_med", "dp"]).filter(pl.col("dp").is_null() | (pl.col("dp") != 75)).collect()
        cods = d["cb"].unique().to_list()
        gt = {c: (fb(c) or 0) if (c and c[0].isalpha()) else 0 for c in cods}
        cie10 = pl.col("cb").str.slice(0, 1).str.contains("[A-Z]").fill_null(False)
        c3 = pl.col("cb").str.slice(0, 3)
        ext = (cie10 & (c3 >= "V01") & (c3 <= "Y98")).fill_null(False)
        d = d.with_columns(pl.col("cb").replace_strict(gt, default=0, return_dtype=pl.Int64).alias("gc"))
        d = d.select(
            "mun", pl.lit(1).alias("n"), ext.cast(pl.Int64).alias("externas"), cie10.cast(pl.Int64).alias("cie10"),
            (pl.col("gc") > 0).cast(pl.Int64).alias("gc_1a4"),
            (pl.col("sit_defun") == 3).fill_null(False).cast(pl.Int64).alias("dom"),
            (pl.col("asis_med") == 2).fill_null(False).cast(pl.Int64).alias("asis_no"),
            ((pl.col("asis_med") == 2) & ext).fill_null(False).cast(pl.Int64).alias("asis_no_ext"))
        partes.append(d.filter(pl.col("mun") > 0).group_by("mun").sum())
    q = pl.concat(partes).with_columns(mun_fusion(pl.col("mun")).alias("mun")).group_by("mun").sum()
    return q.select("mun", pl.col("n").alias("def_total"),
                    (100 * (pl.col("asis_no") - pl.col("asis_no_ext")) / (pl.col("n") - pl.col("externas"))).alias("sin_asistencia"),
                    (100 * pl.col("gc_1a4") / pl.col("cie10")).alias("basura"),
                    (100 * pl.col("dom") / pl.col("n")).alias("casa")).sort("mun")


def basura_sin_hipertension(muns: list[int], anios=ANIOS) -> pl.DataFrame:
    """Proporción de códigos basura entre las defunciones con causa CIE-10 que no son enfermedades hipertensivas
    (I10–I15). Quita el vínculo mecánico con la mortalidad hipertensiva (I10 e I15 son códigos basura de tipo 1)."""
    fb = _clasificador_basura()
    partes = []
    for a in anios:
        d = (_def(a, ["mun", "cb"]).filter(pl.col("mun") > 0).collect()
             .with_columns(mun_fusion(pl.col("mun")).alias("mun"), pl.col("cb").fill_null(""))
             .filter(pl.col("cb").str.slice(0, 1).str.contains("[A-Z]")))
        gt = {c: int((fb(c) or 0) > 0) for c in d["cb"].unique().to_list()}
        d = d.with_columns(pl.col("cb").replace_strict(gt, default=0, return_dtype=pl.Int64).alias("gb"),
                           pl.col("cb").str.slice(0, 3).is_between(pl.lit("I10"), pl.lit("I15")).alias("hta"))
        partes.append(d.group_by("mun").agg(pl.len().alias("n"), pl.col("gb").sum().alias("gb"),
                                            (~pl.col("hta")).sum().alias("n_sin"), pl.col("gb").filter(~pl.col("hta")).sum().alias("gb_sin")))
    t = pl.concat(partes).group_by("mun").sum()
    t = t.with_columns((100 * pl.col("gb") / pl.col("n")).alias("basura_rec"), (100 * pl.col("gb_sin") / pl.col("n_sin")).alias("basura_sin_hta"))
    return pl.DataFrame({"mun": muns}).join(t.select("mun", "basura_rec", "basura_sin_hta"), on="mun", how="left", maintain_order="left")


def completitud_departamental() -> pl.DataFrame:
    """Completitud del registro de defunciones 2015–2017 por departamento: registradas / esperadas, con
    esperadas = Σ tasa de mortalidad del DANE (edad simple, sexo, departamento; base CNPV 2018) × población.
    Misma regla que `scripts/35_completitud_cnpv2018.py` del repositorio EEVV (de donde la toma el estudio 55):
    las defunciones de edad desconocida se reparten en proporción y se excluyen las de sexo desconocido."""
    from python_calamine import CalamineWorkbook
    r = CalamineWorkbook.from_path(str(rutas.RAW / "poblacion" / "dane_mor" / "DCD-Mor-EstSexNal-Dep-1985-2017_VP.xlsx")).get_sheet_by_name("Mx").to_python()
    cab = r[6]; edades = [int(float(x)) for x in cab[5:]]
    filas = []
    for x in r[7:]:
        if x[0] == "" or x[3] != "TOTAL":
            continue
        dp = 0 if x[1] == "00" else int(x[1])
        if not 2015 <= int(x[0]) <= 2017 or dp == 0:
            continue
        for e, m in zip(edades, x[5:]):
            filas.append((int(x[0]), dp, 1 if x[4] == "Hombres" else 2, e, float(m)))
    mx = pl.DataFrame(filas, schema=["anio", "dp", "sexo", "edad", "mx"], orient="row")
    pob = (pl.read_parquet(rutas.PROCESADO / "poblacion_eevv.parquet").filter(pl.col("anio").is_between(2015, 2017) & (pl.col("dp") > 0))
           .group_by(["anio", "dp", "sexo", "edad"]).agg(pl.col("pob").sum())
           .with_columns(pl.col("anio").cast(pl.Int64), pl.col("dp").cast(pl.Int64), pl.col("sexo").cast(pl.Int64), pl.col("edad").cast(pl.Int64)))
    esp = suma_ordenada(mx.join(pob, on=["anio", "dp", "sexo", "edad"], how="inner"), "dp", pl.col("mx") * pl.col("pob"), "E")
    base = (pl.scan_parquet(rutas.PROCESADO / "defunciones_eevv.parquet").select(["anio", "dp", "sexo"])
            .filter(pl.col("anio").is_between(2015, 2017) & pl.col("dp").is_not_null() & (pl.col("dp") != 75) & pl.col("sexo").is_in([1, 2]))
            .group_by("dp").agg(pl.len().alias("D")).collect().with_columns(pl.col("dp").cast(pl.Int64)))
    t = esp.join(base, on="dp", how="left", maintain_order="left").with_columns(pl.col("D").fill_null(0))
    return t.select("dp", (100 * (pl.col("D") / pl.col("E")).round(3)).alias("completitud")).sort("dp")


def departamental_2022() -> pl.DataFrame:
    """Observados y esperados de 2022 por departamento de residencia para las 14 series (estandarización indirecta con
    los mismos 19 grupos de edad), para comparar con la I de Moran departamental del trabajo previo."""
    d = (_def(2022, ["dp", "g", "c667"]).filter(pl.col("dp").is_not_null() & (pl.col("g") >= 0) & ~pl.col("dp").is_in([75, 88]))
         .with_columns(pl.min_horizontal(pl.col("g"), pl.lit(18)).cast(pl.Int8).alias("g")).collect())
    pob = (pl.read_parquet(rutas.PROCESADO / "poblacion_eevv.parquet").filter((pl.col("anio") == 2022) & (pl.col("dp") > 0) & (pl.col("dp") != 88))
           .with_columns(grupo_edad(pl.col("edad")).alias("g")).group_by(["dp", "g"]).agg(pl.col("pob").sum()))
    filas = []
    for k in CLAVES:
        sub = d if k == "todas" else d.filter(pl.col("c667") == k)
        o = sub.group_by(["dp", "g"]).agg(pl.len().alias("o"))
        lam = o.group_by("g").agg(pl.col("o").sum()).join(pob.group_by("g").agg(pl.col("pob").sum()), on="g", how="right") \
               .with_columns((pl.col("o").fill_null(0) / pl.col("pob")).alias("lam")).select("g", "lam")
        e = suma_ordenada(pob.join(lam, on="g"), "dp", pl.col("lam") * pl.col("pob"), "E")
        ob = o.group_by("dp").agg(pl.col("o").sum().alias("O"))
        filas.append(e.join(ob, on="dp", how="left", maintain_order="left").with_columns(pl.col("O").fill_null(0), pl.lit(clave(k)).alias("causa")))
    return pl.concat(filas).select(pl.col("dp").cast(pl.Int32), "causa", "O", "E").sort(["causa", "dp"])


def hav(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371.0088 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def covariables(muns, rural, pob, calidad, completitud, basura_hta) -> pl.DataFrame:
    """Covariables municipales (tabla 1 del manuscrito) de fuentes abiertas: CNPV 2018 integrado al MGN 2018, NBI,
    IPM y etnia 2018 (DANE), DIVIPOLA e IDEAM (estudio 52), REPS 2022, PDET, calidad del registro y completitud."""
    g = json.loads((rutas.EXT63 / "mgn2018_mpio_cnpv.geojson").read_text())
    filas = []
    for f in g["features"]:
        p = f["properties"]
        conoc = sum(p[k] or 0 for k in ("STP51_PRIM", "STP51_SECU", "STP51_SUPE", "STP51_POST", "STP51_13_E"))
        inet = (p["STP19_INT1"] or 0) + (p["STP19_INT2"] or 0)
        filas.append({"mun": int(p["MPIO_CDPMP"]), "area_km2": p["AREA"] / 1e6,
                      "internet": 100 * (p["STP19_INT1"] or 0) / inet if inet else None,
                      "superior": 100 * ((p["STP51_SUPE"] or 0) + (p["STP51_POST"] or 0)) / conoc if conoc else None,
                      "personas_cnpv": p["STP27_PERS"], "lat_c": p["LATITUD"], "lon_c": p["LONGITUD"]})
    mgn = pl.DataFrame(filas)
    nbi = pl.DataFrame([{"mun": int(x["attributes"]["MPIO_CCDGO"]), "nbi": x["attributes"]["NBIC_Total_Prop_de_Personas_en_NBIPorc"]}
                        for x in json.loads((rutas.EXT63 / "nbi2018_mpio.json").read_text())["features"]])
    ipm = pl.DataFrame([{"mun": int(x["attributes"]["MPIO_CCDGO"]), "ipm": x["attributes"]["IPM"]}
                        for x in json.loads((rutas.EXT63 / "ipm2018_mpio.json").read_text())["features"]])
    # celdas vacías en el servicio del DANE = ninguna persona del grupo (también viene vacío el total de unidades)
    et = pl.DataFrame([{"mun": int(x["attributes"]["U_MPIO"]),
                        "indigena": x["attributes"]["PARTICIPACION_PORCENTUAL_PA1_GRP_ETNIC_1"],
                        "afro": x["attributes"]["PARTICIPACION_PORCENTUAL_PA1_GRP_ETNIC_5"]}
                       for x in json.loads((rutas.EXT63 / "etnia2018_mpio.json").read_text())["features"]])
    dv = pl.read_csv(rutas.EXT52 / "divipola_gdxc-w37w.csv", infer_schema=False)
    dv.columns = ["cdp", "ndp", "cod", "nmun", "tipo", "lon", "lat"]
    dv = dv.select(pl.col("cod").cast(pl.Int32).alias("mun"), pl.col("nmun"), pl.col("ndp"),
                   pl.col("lat").str.replace(",", ".").cast(pl.Float64), pl.col("lon").str.replace(",", ".").cast(pl.Float64))
    est = pl.read_csv(rutas.EXT52 / "ideam_estaciones.csv", infer_schema=False).select(
        pl.col("altitud").str.replace(",", ".").cast(pl.Float64, strict=False).alias("alt"),
        pl.col("latitud").str.replace(",", ".").cast(pl.Float64, strict=False).alias("la"),
        pl.col("longitud").str.replace(",", ".").cast(pl.Float64, strict=False).alias("lo")).drop_nulls()
    est = est.filter(pl.col("la").is_between(-5, 14) & pl.col("lo").is_between(-82, -66) & pl.col("alt").is_between(0, 5500))
    ela, elo, ealt = est["la"].to_numpy(), est["lo"].to_numpy(), est["alt"].to_numpy()
    # punto de referencia: DIVIPOLA; si falta (94663 Mapiripana no está en la DIVIPOLA de 2024), el centroide del MGN 2018
    dv = (pl.DataFrame({"mun": muns}).join(dv, on="mun", how="left", maintain_order="left").join(mgn.select("mun", "lat_c", "lon_c"), on="mun", how="left", maintain_order="left")
          .with_columns(pl.col("lat").fill_null(pl.col("lat_c")), pl.col("lon").fill_null(pl.col("lon_c")),
                        pl.col("lat").is_null().alias("coord_mgn")).drop("lat_c", "lon_c"))
    lat, lon = dv["lat"].to_numpy(), dv["lon"].to_numpy()
    alt, kmest = [], []
    for i in range(dv.height):
        d = hav(lat[i], lon[i], ela, elo); j = int(d.argmin()); alt.append(float(ealt[j])); kmest.append(float(d[j]))
    dv = dv.with_columns(pl.Series("altitud", alt), pl.Series("km_estacion", kmest))
    rp = pl.read_csv(rutas.EXT63 / "reps_capacidad_2022.csv", infer_schema=False).with_columns(
        pl.col("c_digo_sede").str.zfill(10).str.slice(0, 5).cast(pl.Int32).alias("mun"),
        pl.col("num_cantidad_capacidad_instalada").cast(pl.Float64, strict=False).fill_null(0).alias("q"))
    rp = rp.with_columns(mun_fusion(pl.col("mun")).alias("mun"))
    uci_desc = pl.col("nom_descripcion_capacidad").str.contains("(?i)intensiv")
    reps = rp.group_by("mun").agg(
        pl.col("q").filter(pl.col("nom_grupo_capacidad") == "CAMAS").sum().alias("camas"),
        pl.col("q").filter((pl.col("nom_grupo_capacidad") == "CAMAS") & uci_desc).sum().alias("camas_uci"))
    cl = pl.read_csv(rutas.CONFIG / "clasificacion_municipal.csv", infer_schema=False).select(
        pl.col("cod_mpio").cast(pl.Int32).alias("mun"), (pl.col("pdet") == "si").cast(pl.Int8).alias("pdet"))
    pobt = pob.group_by("mun").agg((pl.col("pob").sum() / len(ANIOS)).alias("pob_media"))
    b = (pl.DataFrame({"mun": muns}).with_columns((pl.col("mun") // 1000).alias("dp"))
         .join(dv, on="mun", how="left", maintain_order="left").join(mgn, on="mun", how="left", maintain_order="left").join(nbi, on="mun", how="left", maintain_order="left")
         .join(ipm, on="mun", how="left", maintain_order="left").join(et, on="mun", how="left", maintain_order="left").join(reps, on="mun", how="left", maintain_order="left")
         .join(cl, on="mun", how="left", maintain_order="left").join(calidad, on="mun", how="left", maintain_order="left")
         .join(completitud.with_columns(pl.col("dp").cast(pl.Int64)), on="dp", how="left", maintain_order="left")
         .join(rural, on="mun", how="left", maintain_order="left").join(pobt, on="mun", how="left", maintain_order="left")
         .join(basura_hta.select("mun", "basura_sin_hta"), on="mun", how="left", maintain_order="left"))
    b = b.with_columns(pl.col("camas").fill_null(0), pl.col("camas_uci").fill_null(0), pl.col("pdet").fill_null(0),
                       pl.col("indigena").fill_null(0.0), pl.col("afro").fill_null(0.0), (100 * pl.col("rural")).alias("rural"))
    la, lo = b["lat"].to_numpy(), b["lon"].to_numpy()
    cod = b["mun"].to_numpy(); dp = b["dp"].to_numpy()
    idx = {c: i for i, c in enumerate(cod)}
    uci = np.where(b["camas_uci"].to_numpy() > 0)[0]
    dcap, duci = [], []
    for i in range(len(cod)):
        jc = idx.get(dp[i] * 1000 + 1)
        dcap.append(0.0 if jc is None or jc == i else float(hav(la[i], lo[i], la[jc], lo[jc])))
        duci.append(float(hav(la[i], lo[i], la[uci], lo[uci]).min()))
    b = b.with_columns(pl.Series("dist_capital", dcap), pl.Series("dist_uci", duci))
    b = b.with_columns(
        (pl.col("pob_media") / pl.col("area_km2")).log().alias("ldens"),
        (pl.col("dist_capital") + 1).log().alias("ldist_capital"),
        (pl.col("dist_uci") + 1).log().alias("ldist_uci"),
        (1000 * pl.col("camas") / pl.col("pob_media") + 1).log().alias("lcamas"),
        (pl.col("dp").replace_strict(MACRO, default="Otra")).alias("macro"))
    return b.sort("mun")


def redondear(x: float, k: int = 9) -> float:
    return float(round(x, k)) if x is not None and not math.isnan(x) else x
