#!/usr/bin/env python3
"""Verifica cada cifra del manuscrito A3 (manuscrito/manuscrito.md) contra las salidas del repositorio: output/*.json,
output/tablas/*.csv, data/intermedios/ y config/ (cifras citadas de otros trabajos).
Redondea como el texto (mitad hacia arriba) y comprueba que el fragmento citado aparezca tal cual en
manuscrito.md. Termina en «Todo coincide.» o lista las fallas.

    python pruebas/verificar_cifras.py
"""
from __future__ import annotations

import json
import math
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
import polars as pl

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
from agrupamiento import figuras, pesos  # noqa: E402
SAL = RAIZ / "output"
INT = RAIZ / "data" / "intermedios"
_RUTA_TXT = RAIZ / "manuscrito" / "manuscrito.md"
if not _RUTA_TXT.exists():
    print("El texto del artículo no está en este repositorio; se agrega cuando se publique. Las pruebas de cálculo (pytest) sí corren.")
    sys.exit(0)
TXT = (RAIZ / "manuscrito" / "manuscrito.md").read_text(encoding="utf-8")
R = json.loads((SAL / "resultados.json").read_text()); C = R["causas"]
RV = json.loads((SAL / "revision.json").read_text())
BY = json.loads((SAL / "bym2.json").read_text())["causas"]
PA = json.loads((SAL / "panel.json").read_text())
T1 = pl.read_csv(SAL / "tablas" / "tabla1_covariables.csv")
T2 = pl.read_csv(SAL / "tablas" / "tabla2_moran_modelos.csv")
T4 = pl.read_csv(SAL / "tablas" / "tabla4_robustez_W.csv")
T6 = pl.read_csv(SAL / "tablas" / "tabla6_mco_conley.csv")
T8 = pl.read_csv(SAL / "tablas" / "tabla8_efectos_coprincipales.csv")
META = json.loads((INT / "datos_meta.json").read_text())
COV = pl.read_parquet(INT / "covariables.parquet").sort("mun")
OE = pl.read_parquet(INT / "observados_esperados.parquet")
CIT = {r["clave"]: r["valor_texto"] for r in pl.read_csv(RAIZ / "config" / "cifras_citadas.csv").iter_rows(named=True)}
SEL = pl.read_csv(RAIZ / "config" / "seleccion_causas_departamental.csv")
DEP = json.loads((SAL / "departamental.json").read_text())["causas"]
COVS = R["meta"]["covariables"]
fallas, n_ok = [], 0


def es(x, k=2):
    s = str(Decimal(str(float(x))).quantize(Decimal(1).scaleb(-k), rounding=ROUND_HALF_UP))
    neg = s.startswith("-"); s = s.lstrip("-")
    e, _, d = s.partition(".")
    out = f"{int(e):,}".replace(",", ".") + ("," + d if d else "")
    if neg and float(out.replace(".", "").replace(",", ".")) != 0:
        out = "−" + out
    return out


def chk(que, valor, texto, k=2, contexto=None):
    """Compara el valor recalculado (redondeado a k decimales) con la cifra del texto y busca el contexto."""
    global n_ok
    calc = valor if isinstance(valor, str) else es(valor, k)
    ctx = texto if contexto is None else contexto
    if calc == texto and ctx in TXT:
        n_ok += 1
    else:
        fallas.append(f"{que}: calculado {calc}, texto {texto}, contexto {'encontrado' if ctx in TXT else 'NO encontrado'}: «{ctx[:90]}»")


def si(cond):
    return "sí" if cond else "no"


def fila(t, **kw):
    for k, v in kw.items():
        t = t.filter(pl.col(k) == v)
    return t.row(0, named=True)


pct = lambda x: 100 * math.expm1(x)
DOCE = [c for c in C if C[c]["moran_y"]["I"] > 0.1]

# ── datos ────────────────────────────────────────────────────────────────
p = META["perdidas"]
tot = sum(v["total"] for v in p.values()); sinmun = sum(v["sin_mun"] for v in p.values())
sa = sum(v["san_andres"] for v in p.values()); se_ = sum(v["sin_edad"] for v in p.values())
an = int(fila(T2, causa="todas")["defunciones"])
chk("defunciones", tot, "2.260.923", 0, "De 2.260.923 defunciones del periodo, 26.407 (1,2 %) no tienen municipio de residencia, 2.659 (0,12 %) son de San Andrés y Providencia y 3 no tienen edad; quedan 2.231.854")
chk("sin municipio", sinmun, "26.407", 0); chk("sin municipio %", 100 * sinmun / tot, "1,2", 1)
chk("San Andrés", sa, "2.659", 0); chk("San Andrés %", 100 * sa / tot, "0,12", 2, "Como pesan 0,12 % de las defunciones")
chk("sin edad", se_, "3", 0); chk("analizadas", an, "2.231.854", 0); chk("cuadre", tot - sinmun - sa - se_, "2.231.854", 0)
chk("millones", an / 1e6, "2,23", 2, "(2,23 millones de defunciones, 2018–2025, 1.120 municipios)")
chk("municipios", R["meta"]["n"], "1.120", 0)
chk("ceros vectores", C["103"]["ceros"], "556", 0, "556 de los 1.120 no registraron")
chk("vecinos", R["meta"]["W"]["reina"]["vecinos_media"], "5,7", 1, "(5,7 vecinos en promedio; mínimo 1, máximo 19)")
chk("vecinos min", R["meta"]["W"]["reina"]["vecinos_min"], "1", 0); chk("vecinos max", R["meta"]["W"]["reina"]["vecinos_max"], "19", 0)
chk("forzados", R["meta"]["W"]["invdist"]["con_vecino_forzado"], "11", 0, "los 11 municipios sin vecinos a esa distancia")
chk("UCI", int((COV["camas_uci"] > 0).sum()), "93", 0, "93 municipios las tienen")
chk("PDET", int(COV["pdet"].sum()), "170", 0, "PDET, 170 municipios")
chk("km estación", float(COV["km_estacion"].median()), "0,8", 1, "a una mediana de 0,8 km")
W, _ = pesos.todas(COV); Wd = W["reina"].full()[0]
z = lambda v: (v - v.mean()) / v.std()
corrlag = {v: float(np.corrcoef(z(COV[v].to_numpy().astype(float)), Wd @ z(COV[v].to_numpy().astype(float)))[0, 1]) for v in COVS}
chk("completitud valores", int(COV["completitud"].n_unique()), "32", 0, "La completitud tiene solo 32 valores distintos y una correlación de 0,94 con su propio rezago espacial")
chk("corr rezago completitud", corrlag["completitud"], "0,94", 2, "va de 0,17 en las camas a 0,94 en la completitud")
chk("corr rezago camas", corrlag["lcamas"], "0,17", 2)
chk("corr rezago extremos", f"{min(corrlag, key=corrlag.get)}/{max(corrlag, key=corrlag.get)}", "lcamas/completitud", 0, "va de 0,17 en las camas a 0,94 en la completitud")
Iv = RV["I_covariables"]
chk("I cov camas", Iv["lcamas"], "0,09", 2, "las I reales van de 0,09 en las camas a 0,85 en la completitud")
chk("I cov completitud", Iv["completitud"], "0,85", 2)
chk("I cov extremos", f"{min(Iv, key=Iv.get)}/{max(Iv, key=Iv.get)}", "lcamas/completitud", 0, "las I reales van de 0,09 en las camas a 0,85 en la completitud")


def contr(c):
    d = COV.select("mun").join(OE.filter(pl.col("causa") == c), on="mun", how="left")
    O, E = d["O"].to_numpy().astype(float), d["E"].to_numpy()
    r = O / E; m = O.sum() / E.sum(); s2 = max((E * (r - m) ** 2).sum() / E.sum() - m / E.mean(), 0)
    return s2 / (s2 + m / E)


chk("C ~ ruralidad estómago", float(np.corrcoef(contr("201"), COV["rural"].to_numpy())[0, 1]), "−0,52", 2,
    "(correlación de −0,52 en el cáncer de estómago) o la altitud (−0,36 en la dificultad respiratoria del RN)")
chk("C ~ altitud RN", float(np.corrcoef(contr("404"), COV["altitud"].to_numpy())[0, 1]), "−0,36", 2)

SE = figuras.smr_eb().join(COV.select("mun", "macro"), on="mun")


def reg(c, m):
    q = SE.filter((pl.col("causa") == c) & (pl.col("macro") == m))
    return float(q["O"].sum() / q["E"].sum())


chk("estómago Caribe", reg("201", "Caribe"), "0,49", 2, "en el Caribe (razón de 0,49) y lo superó en la región andina (1,12)")
chk("estómago Andina", reg("201", "Andina"), "1,12", 2)
chk("homicidios Pacífica", 100 * (reg("512", "Pacífica") - 1), "79", 0, "los homicidios superaron en 79 % el promedio en la región Pacífica")
chk("vectores Orinoquía", reg("103", "Orinoquía"), "4,28", 2, "más de cuatro veces el promedio en la Orinoquía (4,28)")
cd = COV.group_by("dp").agg(pl.col("completitud").first()).sort("completitud")
chk("completitud Vichada", float(cd["completitud"][0]), "29,3", 1, "va de 29,3 % en Vichada a 89,3 % en Bogotá")
chk("completitud Bogotá", float(cd["completitud"][-1]), "89,3", 1)
chk("depto mín y máx", f"{cd['dp'][0]}/{cd['dp'][-1]}", "99/11", 0, "va de 29,3 % en Vichada a 89,3 % en Bogotá")
chk("cifras del estudio 55", "/".join(CIT[k] for k in ("certifica_no_tratante_2023_2025", "sin_asistencia_enfermedad_2023_2025", "sin_asistencia_narino", "codigos_basura_2023_2025")),
    "54,8/17,0/32,9/11,3", 0,
    "el 54,8 % de las muertes las certificó un médico distinto del tratante, que 17,0 % de las muertes por enfermedad ocurrió sin asistencia médica (32,9 % en Nariño)")

# tabla 1
for v, nom, vals in [("nbi", "NBI 2018 (%)", ["22,9", "17,8", "7,4", "16,8", "47,0"]), ("rural", "Población rural (%)", ["55,5", "23,7", "19,5", "58,9", "84,1"]),
                     ("indigena", "Población indígena (%)", ["7,4", "19,2", "0,0", "0,1", "21,8"]), ("afro", "Población afrocolombiana (%)", ["6,9", "19,1", "0,1", "0,4", "19,0"]),
                     ("internet", "Viviendas con internet (%)", ["12,2", "14,3", "1,9", "6,9", "28,3"]), ("superior", "Educación superior (%)", ["10,1", "9,0", "2,8", "7,1", "22,3"]),
                     ("sin_asistencia", "Muertes por enfermedad sin asistencia (%)", ["27,6", "11,4", "13,9", "26,3", "42,7"]),
                     ("basura", "Códigos basura (%)", ["12,7", "5,3", "8,0", "11,8", "18,2"]), ("casa", "Muertes en casa (%)", ["31,8", "10,1", "19,5", "31,3", "45,0"]),
                     ("completitud", "Completitud departamental DANE (%)", ["73,2", "11,6", "61,1", "75,6", "85,4"])]:
    r = fila(T1, variable=v)
    chk(f"tabla 1 {v}", " | ".join(es(r[c], 1) for c in ("media", "de", "p10", "mediana", "p90")), " | ".join(vals), 0,
        f"| {nom} | {r['bloque']} | " + " | ".join(vals) + " |")
r = fila(T1, variable="altitud")
chk("tabla 1 altitud", " | ".join(es(r[c], 0) for c in ("media", "de", "p10", "mediana", "p90")), "1.119 | 938 | 25 | 1.046 | 2.549", 0, "| Altitud (m) | G | 1.119 | 938 | 25 | 1.046 | 2.549 |")
chk("tabla 1 UCI", float(COV["dist_uci"].median()), "27,3", 1, "| Distancia a UCI (km, mediana) | G | | | | 27,3 | |")
chk("tabla 1 PDET", 100 * float(COV["pdet"].mean()), "15,2", 1, "| Municipio PDET (%) | G | 15,2 |")
chk("DE altitud en el texto", r["de"], "938", 0, "Una desviación estándar más de altitud (938 m)")

# ── autocorrelación, placebo y tabla 2 ───────────────────────────────────
chk("E[I] SGQ", RV["EI_H0"]["SGQ"], "−0,007", 3, "−0,007 con los tres bloques y −0,025 si se añaden efectos fijos de departamento")
chk("E[I] SGQ+FE", RV["EI_H0"]["SGQ+FE"], "−0,025", 3)
chk("E[I] en tabla 2", RV["EI_H0"]["SGQ+FE"], "−0,025", 3, "el valor esperado de esta I bajo independencia es −0,025")
chk("E[I] en 5.2", RV["EI_H0"]["SGQ+FE"], "−0,025", 3, "frente a un valor esperado de −0,025 bajo independencia")
chk("pasos no monótonos", sum(len(d["no_monotono"]) for d in RV["causas"].values()), "15", 0, "en 15 de los 168 pasos posibles (14 series por 12 pasos)")
chk("causas no monótonas", ",".join(sorted(c for c, d in RV["causas"].items() if d["no_monotono"])), "103,108,204,210,406,613", 0,
    "Ocurre en las malformaciones, las otras afecciones perinatales, el cáncer de pulmón, el de próstata, la neumonía y las enfermedades por vectores")
chk("réplicas", RV["meta"]["replicas"], "200", 0, "Con 200 réplicas se calculan")
NOMB = {"512": "Homicidios", "201": "Cáncer de estómago", "103": "Enf. por vectores y rabia", "204": "Cáncer de pulmón", "108": "Neumonía e IRA",
        "302": "Enf. hipertensivas", "todas": "Todas las causas", "505": "Ahogamiento", "604": "Alzheimer, Parkinson y otras", "210": "Cáncer de próstata",
        "404": "Dificultad respiratoria del RN", "212": "Leucemia", "406": "Otras afecciones perinatales", "613": "Malformaciones congénitas"}
for c, nom in NOMB.items():
    d = RV["causas"][c]; r2 = fila(T2, causa=c)
    cel = [es(r2["defunciones"], 0), es(d["I0"], 2), es(100 * d["caida_real"], 1), es(100 * d["caida_placebo_media"], 1), es(100 * d["caida_exceso"], 1),
           es(d["p_placebo"], 3), es(d["I_SGQFE"], 2), es(r2["RLMe"], 1), es(r2["RLMl"], 1)]
    linea = f"| {nom} | " + " | ".join(cel) + " |"
    chk(f"tabla 2 {c}", linea, linea, 0, linea)
    chk(f"I inicial igual en las dos salidas {c}", si(abs(d["I0"] - C[c]["moran_y"]["I"]) < 1e-9), "sí", 0)
chk("p permutación", max(C[c]["moran_y"]["p_sim"] for c in C), "0,001", 3,
    "La I inicial tiene p = 0,001 con 999 permutaciones en todas las series.")
chk("Assunção-Reis", f"{es(C['512']['moran_ar']['I'])}, {es(C['201']['moran_ar']['I'])} y {es(C['todas']['moran_ar']['I'])}", "0,59, 0,51 y 0,29", 0,
    "(0,59, 0,51 y 0,29 en homicidios, estómago y total)")
chk("I inicial textual", f"{es(C['512']['moran_y']['I'])}/{es(C['201']['moran_y']['I'])}/{es(C['103']['moran_y']['I'])}/{es(C['204']['moran_y']['I'])}/{es(C['todas']['moran_y']['I'])}",
    "0,62/0,58/0,48/0,44/0,29", 0, "Los homicidios (I = 0,62), el cáncer de estómago (0,58), las enfermedades por vectores (0,48) y el cáncer de pulmón (0,44) son las series más agrupadas; la mortalidad total tiene I = 0,29")
chk("12 series", len(DOCE), "12", 0, "en las 12 series con I inicial mayor de 0,1")
chk("rango I 12 series", f"{es(min(C[c]['moran_y']['I'] for c in DOCE))} a {es(max(C[c]['moran_y']['I'] for c in DOCE))}", "0,12 a 0,62", 0, "(I de 0,12 a 0,62)")
cr = {c: RV["causas"][c]["caida_real"] for c in DOCE}
chk("caída real", f"{es(100 * min(cr.values()), 0)} % ({min(cr, key=cr.get)}) y {es(100 * max(cr.values()), 0)} % ({max(cr, key=cr.get)})", "21 % (302) y 55 % (404)", 0,
    "Las covariables reales reducen la I entre 21 % (enfermedades hipertensivas) y 55 % (dificultad respiratoria del RN)")
chk("caída real en la introducción", f"{es(100 * min(cr.values()), 0)} % y {es(100 * max(cr.values()), 0)} %", "21 % y 55 %", 0, "las covariables reales la reducen entre 21 % y 55 %")
cp = {c: RV["causas"][c]["caida_placebo_media"] for c in DOCE}
chk("placebo", f"{es(100 * min(cp.values()), 0)} % ({min(cp, key=cp.get)}) y {es(100 * max(cp.values()), 0)} % ({max(cp, key=cp.get)})", "10 % (todas) y 22 % (404)", 0,
    "reducen en promedio entre 10 % (mortalidad total) y 22 % (dificultad respiratoria del RN)")
chk("placebo en el resumen", f"{es(100 * min(cp.values()), 0)} % y {es(100 * max(cp.values()), 0)} %", "10 % y 22 %", 0, "reducen la I de Moran entre 10 % y 22 %")
chk("placebo en la introducción", f"{es(100 * min(cp.values()), 0)} % y {es(100 * max(cp.values()), 0)} %", "10 % y 22 %", 0, "16 covariables placebo la reducen entre 10 % y 22 %")
ok11 = [c for c in DOCE if RV["causas"][c]["p_placebo"] <= 0.03]
chk("series con exceso", len(ok11), "11", 0, "en 11 de esas 12 series, con p de 0,03 o menos")
chk("la excluida es hipertensivas", ",".join(c for c in DOCE if c not in ok11), "302", 0)
ex = {c: RV["causas"][c]["caida_exceso"] for c in ok11}
chk("exceso", f"{es(100 * min(ex.values()), 1)} puntos ({min(ex, key=ex.get)}) a {es(100 * max(ex.values()), 1)} ({max(ex, key=ex.get)})", "12,5 puntos (210) a 32,8 (404)", 0,
    "El exceso sobre el placebo va de 12,5 puntos (cáncer de próstata) a 32,8 (dificultad respiratoria del RN)")
chk("exceso redondeado", f"{es(100 * min(ex.values()), 0)} y {es(100 * max(ex.values()), 0)}", "13 y 33", 0, "entre 13 y 33 puntos más en 11 de 12 series")
chk("exceso redondeado intro", f"{es(100 * min(ex.values()), 0)} a {es(100 * max(ex.values()), 0)}", "13 a 33", 0, "va de 13 a 33 puntos porcentuales en 11 series")
d302 = RV["causas"]["302"]
chk("exceso hipertensivas", f"{es(100 * d302['caida_exceso'], 1)} puntos (p = {es(d302['p_placebo'], 2)})", "5,9 puntos (p = 0,14)", 0, "En las enfermedades hipertensivas el exceso es 5,9 puntos (p = 0,14)")
for c, nums in (("todas", "63/42/37"), ("512", "59/35/38")):
    d = RV["causas"][c]
    chk(f"numerador {c}", f"{es(100 * d['caida_num'], 0)}/{es(100 * d['caida_den'], 0)}/{es(100 * d['caida_real'], 0)}", nums, 0,
        "el numerador de la I baja 63 % con los tres bloques, el denominador 42 % y la I 37 %; en los homicidios, 59 %, 35 % y 38 %")
shx = lambda c, b: 100 * RV["causas"][c]["bloques"][b]["shapley_exceso"]
shi = lambda c, b: [100 * x for x in RV["causas"][c]["bloques"][b]["shapley_exceso_ic"]]
qpos = sorted(c for c in RV["causas"] if RV["causas"][c]["bloques"]["Q"]["shapley_exceso_ic"][0] > 0)
chk("Q supera al placebo", f"{len(qpos)}: {','.join(qpos)}", "6: 103,201,212,404,512,todas", 0, "El bloque de calidad del registro supera al placebo en seis series")
chk("Q en el resumen", len(qpos), "6", 0, "La calidad del registro aporta más que el placebo en seis series")
chk("Q total", f"{es(shx('todas', 'Q'), 1)} puntos de la I inicial; intervalo placebo de {es(shi('todas', 'Q')[0], 1)} a {es(shi('todas', 'Q')[1], 1)}",
    "14,2 puntos de la I inicial; intervalo placebo de 11,0 a 16,1", 0, "la mortalidad total (14,2 puntos de la I inicial; intervalo placebo de 11,0 a 16,1)")
chk("Q resto", "/".join(es(shx(c, "Q"), 1) for c in ("512", "212", "201", "404", "103")), "12,9/11,6/11,0/9,9/8,1", 0,
    "los homicidios (12,9), la leucemia (11,6), el cáncer de estómago (11,0), la dificultad respiratoria del RN (9,9) y las enfermedades por vectores (8,1)")
chk("Q intro", f"{es(shx('todas', 'Q'), 0)} y {es(shx('512', 'Q'), 0)}", "14 y 13", 0, "entre ellas la mortalidad total (14 puntos) y los homicidios (13 puntos)")
chk("Q hipertensivas", f"{es(shx('302', 'Q'), 1)} puntos ({es(shi('302', 'Q')[0], 1)} a {es(shi('302', 'Q')[1], 1)})", "0,9 puntos (−6,2 a 4,5)", 0,
    "En las enfermedades hipertensivas su exceso es 0,9 puntos (−6,2 a 4,5)")
chk("S por causa", "/".join(es(shx(c, "S"), 1) for c in ("204", "404", "108", "505")), "23,5/23,2/20,4/17,8", 0,
    "pesa más en el cáncer de pulmón (23,5), la dificultad respiratoria del RN (23,2), la neumonía (20,4) y el ahogamiento (17,8)")
chk("G estómago", shx("201", "G"), "15,7", 1, "el de geografía y oferta, en el cáncer de estómago (15,7)")
chk("S total", f"{es(shx('todas', 'S'), 1)}; {es(shi('todas', 'S')[0], 1)} a {es(shi('todas', 'S')[1], 1)}", "8,7; 3,5 a 11,8", 0,
    "(14,2; 11,0 a 16,1) es mayor que el del bloque sociodemográfico (8,7; 3,5 a 11,8)")
ife = {c: RV["causas"][c]["I_SGQFE"] for c in DOCE}
chk("I con EF", f"{es(min(ife.values()))} y {es(max(ife.values()))}", "0,01 y 0,21", 0, "la I residual queda entre 0,01 y 0,21 en las 12 series")
chk("I con EF en el resumen", es(max(ife.values())), "0,21", 0, "Con efectos fijos de departamento la I residual")
chk("I con EF estómago y homicidios", f"{es(ife['201'])} y {es(ife['512'])}", "0,12 y 0,21", 0, "(I residual con efectos fijos de 0,12 y 0,21)")
chk("I municipal perinatales", f"{es(C['613']['moran_y']['I'])} y {es(C['406']['moran_y']['I'])}", "0,06 y 0,07", 0,
    "Las malformaciones congénitas y las otras afecciones perinatales tienen I municipal de 0,06 y 0,07.")
I22 = {str(r["codigo_667"]): r["I_moran_2022"] for r in SEL.iter_rows(named=True)}
chk("I departamental del trabajo previo", f"{es(I22['613'])} y {es(I22['406'])}", "0,37 y 0,31", 0, "El trabajo previo les atribuía una I departamental de 0,37 y 0,31 en 2022")
chk("I departamental recalculada", f"{es(DEP['613']['I'])} y {es(DEP['406']['I'])} (p de {es(DEP['613']['p_perm'])} y {es(DEP['406']['p_perm'])})", "−0,09 y −0,11 (p de 0,69 y 0,76)", 0,
    "la I departamental de 2022 es −0,09 y −0,11 (p de 0,69 y 0,76), y la de la dificultad respiratoria del RN baja de 0,27 a 0,08 (p = 0,12); las otras nueve causas conservan p < 0,05")
chk("RN departamental", f"de {es(I22['404'])} a {es(DEP['404']['I'])} (p = {es(DEP['404']['p_perm'])})", "de 0,27 a 0,08 (p = 0,12)", 0)
doce_sel = [str(c) for c in SEL.filter(pl.col("seleccionada") == 1)["codigo_667"]]
no_dep = sorted(c for c in doce_sel if DEP[c]["p_perm"] >= 0.05)
chk("tres no se agrupan entre departamentos", ",".join(no_dep), "404,406,613", 0,
    "tres de las 12 (malformaciones congénitas, otras afecciones perinatales y dificultad respiratoria del recién nacido) dejan de agruparse entre departamentos")
chk("nueve mantienen p < 0,05", len(doce_sel) - len(no_dep), "9", 0, "nueve de las 12 mantienen p < 0,05 con 999 permutaciones y tres no")
chk("selección del trabajo previo", f"{SEL.height}/{int((SEL['p_9999'] < 0.05).sum())}/{int(SEL['seleccionada'].sum())}", "61/27/12", 0,
    "mantienen q < 0,05 después de la corrección de Benjamini y Hochberg sobre 61 causas")
chk("BYM2 φ nulo perinatales", f"{es(100 * BY['613']['nulo']['phi'], 0)} % y {es(100 * BY['406']['nulo']['phi'], 0)} %", "26 % y 50 %", 0,
    "atribuye 26 % y 50 % de la varianza del riesgo a un componente espacialmente estructurado, y el efecto aleatorio estimado tiene una I de 0,33 y 0,54")
chk("BYM2 I del efecto", f"{es(BY['613']['moran_u_nulo'])} y {es(BY['406']['moran_u_nulo'])}", "0,33 y 0,54", 0)
cs = {c: 1 - BY[c]["completo"]["sigma2"] / BY[c]["nulo"]["sigma2"] for c in DOCE}
chk("BYM2 caída de σ²", f"{es(100 * min(cs.values()), 0)} % ({min(cs, key=cs.get)}) y {es(100 * max(cs.values()), 0)} % ({max(cs, key=cs.get)})", "22 % (302) y 56 % (604)", 0,
    "La varianza marginal del efecto BYM2 cae entre 22 % (enfermedades hipertensivas) y 56 % (Alzheimer)")
fr0 = sorted(c for c in BY if BY[c]["nulo"]["phi"] > 0.98); fr1 = sorted(c for c in BY if BY[c]["completo"]["phi"] > 0.98)
chk("frontera", f"{','.join(fr0)} | {','.join(fr1)}", "103,201,512 | 103,204,210", 0,
    "en el modelo sin covariables del cáncer de estómago, las enfermedades por vectores y los homicidios, y en el modelo con covariables del cáncer de pulmón, el de próstata y las enfermedades por vectores")
chk("seis ajustes en la frontera", len(fr0) + len(fr1), "6", 0, "en seis ajustes $\\varphi$ quedó en la frontera")
chk("BYM2 convergencia", si(all(BY[c][m]["conv"] for c in BY for m in ("nulo", "completo"))), "sí", 0)

# ── forma de la dependencia y panel ──────────────────────────────────────
lag = ["512", "201", "204", "103", "108", "404"]
chk("LM-R rezago", f"{es(min(C[c]['ols']['rlm_lag'][0] for c in lag), 1)} a {es(max(C[c]['ols']['rlm_lag'][0] for c in lag), 1)}", "16,7 a 151,3", 0, "(estadísticos de 16,7 a 151,3)")
chk("rezago mayor que error", si(all(C[c]["ols"]["rlm_lag"][0] > C[c]["ols"]["rlm_error"][0] for c in lag)), "sí", 0)
chk("LM-R total", f"{es(C['todas']['ols']['rlm_error'][0], 1)} frente a {es(C['todas']['ols']['rlm_lag'][0], 1)}", "18,2 frente a 2,5", 0,
    "un error espacial en la mortalidad total (18,2 frente a 2,5)")
fe = {c: C[c]["fe_espacial"] for c in C}
chk("ρ con EF", f"{es(fe['512']['rho_sdm'])} (error estándar {es(fe['512']['se_rho_sdm'])}) en los homicidios, {es(fe['103']['rho_sdm'])} en las enfermedades por vectores y {es(fe['201']['rho_sdm'])}",
    "0,44 (error estándar 0,04) en los homicidios, 0,33 en las enfermedades por vectores y 0,28", 0,
    "sigue siendo 0,44 (error estándar 0,04) en los homicidios, 0,33 en las enfermedades por vectores y 0,28 en el cáncer de estómago")
chk("ρ con EF baja en todas", si(all(fe[c]["rho_sdm"] < C[c]["modelos"]["SDM"]["rho"] for c in C)), "sí", 0, "el $\\rho$ del SDM baja en todas las series")
vals = [fe[c][k] for c in ("613", "212", "406", "404") for k in ("rho_sdm", "lam_sdem")]
chk("rango con EF", f"{es(min(vals))} y {es(max(vals))}", "−0,03 y 0,05", 0, "$\\rho$ y $\\lambda$ quedan entre −0,03 y 0,05")
chk("prueba de rezago del estómago", C["201"]["ols"]["rlm_lag"][0], "109,4", 1, "tiene una de las pruebas de rezago más altas (109,4)")
chk("panel n", PA["512"]["anio"]["n"], "7.834", 0, "7.834 observaciones municipio-año")
chk("panel homicidios", f"{es(PA['512']['anio']['delta'], 3)}, error estándar {es(PA['512']['anio']['se_delta'], 3)}", "0,197, error estándar 0,048", 0,
    "($\\delta$ = 0,197, error estándar 0,048), y la asociación se mantiene con efectos fijos de departamento por año (0,086; 0,034)")
chk("panel homicidios dpto-año", f"{es(PA['512']['dpto_anio']['delta'], 3)}; {es(PA['512']['dpto_anio']['se_delta'], 3)}", "0,086; 0,034", 0)
chk("panel estómago", f"{es(PA['201']['anio']['delta'], 3)}; {es(PA['201']['anio']['se_delta'], 3)}", "0,007; 0,026", 0, "No aparece en el cáncer de estómago (0,007; 0,026)")
chk("panel próstata", f"{es(PA['210']['dpto_anio']['delta'], 3)}; {es(PA['210']['dpto_anio']['se_delta'], 3)}", "0,090; 0,036", 0,
    "en el cáncer de próstata (0,090; 0,036, con efectos de departamento por año)")
for c, nom in (("512", "Homicidios"), ("604", "Alzheimer, Parkinson y otras"), ("108", "Neumonía e IRA"), ("210", "Cáncer de próstata"),
               ("302", "Enf. hipertensivas"), ("103", "Enf. por vectores y rabia"), ("todas", "Todas las causas"), ("201", "Cáncer de estómago")):
    a, b = PA[c]["anio"], PA[c]["dpto_anio"]
    linea = f"| {nom} | {es(a['delta'], 3)} ({es(a['se_delta'], 3)}) | {es(b['delta'], 3)} ({es(b['se_delta'], 3)}) |"
    chk(f"tabla 3 {c}", linea, linea, 0, linea)
chk("panel no estimable", si(all(PA["212"][m]["estimable"] is False for m in ("anio", "dpto_anio")) and
       all(PA[c][m]["estimable"] for c in PA if c != "212" for m in ("anio", "dpto_anio"))), "sí", 0,
    "En la leucemia el suavizado anual contrae todas las SMR al promedio en los ocho años y $\\delta$ no es estimable")
chk("panel perinatales", si(all(abs(PA["406"][m]["delta"]) < 1e-4 and PA["406"][m]["se_delta"] < 1e-4 for m in ("anio", "dpto_anio"))),
    "sí", 0, "en las otras afecciones perinatales lo hace en seis de los ocho años y $\\delta$ es prácticamente cero")

# ── efectos coprincipales (tabla 4) ──────────────────────────────────────
chk("q < 0,05", f"{int((T8['mcp_q_cl'] < 0.05).sum())}/{int((T8['bym_q'] < 0.05).sum())}/{int(T8['sobrevive'].sum())}", "52/28/14", 0,
    "52 tienen q < 0,05 en el SLX ponderado, 28 en el Poisson-BYM2 y 14 cumplen el criterio de supervivencia")
chk("sobreviven en el resumen", int(T8["sobrevive"].sum()), "14", 0, "Sobreviven 14 de 224 asociaciones")
chk("sobreviven en la introducción", int(T8["sobrevive"].sum()), "14", 0, "De 224 asociaciones entre causa y covariable, 14 tienen el mismo signo")
chk("combinaciones", T8.height, "224", 0)
chk("q Conley", int((T8["mcp_q_conley"] < 0.05).sum()), "60", 0, "los errores de Conley dan 60 combinaciones con q < 0,05, frente a 52 con errores agrupados")
sobrev = sorted(f"{r['causa']}:{r['variable']}" for r in T8.filter(pl.col("sobrevive")).iter_rows(named=True))
esper = sorted(["todas:rural", "604:completitud", "201:rural", "201:altitud", "201:pdet", "204:rural", "108:rural", "108:indigena", "404:nbi", "404:basura",
                "512:superior", "512:pdet", "512:sin_asistencia", "512:basura"])
chk("lista de sobrevivientes", ",".join(sobrev), ",".join(esper), 0,
    "Las otras ocho asociaciones que sobreviven son de la ruralidad (negativa en la mortalidad total, el cáncer de estómago, el de pulmón y la neumonía), de PDET en el cáncer de estómago, de los códigos basura en la dificultad respiratoria del RN y en los homicidios, y de las muertes sin asistencia en los homicidios")
chk("ruralidad negativa", si(all(r["mcp_total"] < 0 and r["bym_total"] < 0 for r in T8.filter(pl.col("sobrevive") & (pl.col("variable") == "rural")).iter_rows(named=True))), "sí", 0)
F = {(f["causa"], f["variable"]): f for f in figuras.efectos_seleccionados(R, RV, BY)}
ROT = {("201", "altitud"): "Cáncer de estómago · altitud", ("404", "nbi"): "Dificultad resp. del RN · NBI", ("108", "indigena"): "Neumonía e IRA · indígena",
       ("604", "completitud"): "Alzheimer y otras · completitud", ("512", "superior"): "Homicidios · educación superior", ("512", "pdet"): "Homicidios · PDET",
       ("103", "altitud"): "Enf. por vectores · altitud", ("204", "internet"): "Cáncer de pulmón · internet", ("210", "afro"): "Cáncer de próstata · afrocolombiana",
       ("604", "internet"): "Alzheimer y otras · internet", ("todas", "sin_asistencia"): "Todas las causas · sin asistencia",
       ("302", "basura"): "Enf. hipertensivas · códigos basura", ("302", "basura_sin_hta"): "Enf. hipertensivas · basura sin I10–I15"}
for (c, v), rot in ROT.items():
    f = F[(c, v)]; m, b = f["mcp"], f["bym"]
    if v == "basura_sin_hta":
        q, sv = "", "No"
    else:
        q = "< 0,001" if f["mcp_q"] < 0.001 else es(f["mcp_q"], 3)
        sv = "Sí" if fila(T8, causa=c, variable=v)["sobrevive"] else "No"
    linea = f"| {rot} | {es(m[0], 1)} ({es(m[1], 1)} a {es(m[2], 1)}) | {q} | {es(b[0], 1)} ({es(b[1], 1)} a {es(b[2], 1)}) | {sv} |".replace("|  |", "| |")
    chk(f"tabla 4 {c} {v}", linea, linea, 0, linea)
chk("altitud en el texto", f"{es(F[('201', 'altitud')]['mcp'][0], 1)} % más mortalidad por cáncer de estómago en el SLX ponderado y {es(F[('201', 'altitud')]['bym'][0], 1)} %",
    "24,6 % más mortalidad por cáncer de estómago en el SLX ponderado y 28,9 %", 0, "se asocia con 24,6 % más mortalidad por cáncer de estómago en el SLX ponderado y 28,9 % en el modelo de conteo")
chk("Chow altitud", C["201"]["regimenes"]["chow_vars"]["altitud"][1], "0,19", 2, "no rechaza que el coeficiente sea igual en las cinco regiones (p = 0,19)")
fb, fs = F[("302", "basura")], F[("302", "basura_sin_hta")]
chk("códigos basura sin I10–I15", f"{es(fb['mcp'][0], 1)} % a {es(fs['mcp'][0], 1)} % ({es(fs['mcp'][1], 1)} a {es(fs['mcp'][2], 1)}) en el SLX ponderado y de {es(fb['bym'][0], 1)} % a {es(fs['bym'][0], 1)} % ({es(fs['bym'][1], 1)} a {es(fs['bym'][2], 1)})",
    "18,9 % a 7,6 % (−6,5 a 23,8) en el SLX ponderado y de 10,8 % a −2,0 % (−11,9 a 9,0)", 0,
    "el efecto total pasa de 18,9 % a 7,6 % (−6,5 a 23,8) en el SLX ponderado y de 10,8 % a −2,0 % (−11,9 a 9,0) en el modelo de conteo")
chk("sin asociación con el indicador corregido", si(fs["mcp"][1] < 0 < fs["mcp"][2] and fs["bym"][1] < 0 < fs["bym"][2]), "sí", 0, "Con el indicador corregido no queda asociación")
chk("PDET total", F[("512", "pdet")]["mcp"][0], "86,3", 1, "El efecto total de PDET (86,3 % en el SLX ponderado)")
j = COVS.index("pdet"); mdir = RV["slx"]["512"]["MCP"]["dir"][j]; sdir = RV["slx"]["512"]["MCP"]["cl"]["se_dir"][j]
chk("PDET propio", pct(mdir), "−0,2", 1, "el efecto propio del municipio (−0,2 %) no se distingue de cero")
chk("PDET propio incluye cero", si(mdir - 1.96 * sdir < 0 < mdir + 1.96 * sdir), "sí", 0)
chk("educación y homicidios", -F[("512", "superior")]["mcp"][0], "43,5", 1, "La educación superior se asocia con 43,5 % menos homicidios (efecto total)")
rr = fila(T8, causa="404", variable="ldist_uci")
chk("UCI y RN", si(rr["mcp_q_cl"] >= 0.05 and rr["bym_q"] >= 0.05), "sí", 0, "La distancia a cuidado intensivo no se asocia con la dificultad respiratoria del RN en ninguna")
for (c, v) in (("103", "altitud"), ("404", "altitud"), ("204", "internet")):
    r_ = fila(T8, causa=c, variable=v)
    chk(f"una especificación y no la otra {c} {v}", si((r_["mcp_q_cl"] < 0.05) != (r_["bym_q"] < 0.05)), "sí", 0,
        "La altitud y las enfermedades por vectores, la altitud y la dificultad respiratoria del RN y el internet y el cáncer de pulmón aparecen en una especificación pero no en la otra")
for (c, v) in (("210", "afro"), ("604", "internet")):
    chk(f"no sobrevive {c} {v}", si(fila(T8, causa=c, variable=v)["mcp_q_cl"] >= 0.05), "sí", 0,
        "(proporción afrocolombiana y cáncer de próstata, internet y Alzheimer) no tienen efecto total distinto de cero en el SLX ponderado")

# ── robustez y anexo ─────────────────────────────────────────────────────
r = fila(T4, causa="512"); cai = [1 - r[f"I_SGQ_{w}"] / r[f"I_y_{w}"] for w in ("reina", "knn5", "knn8", "invdist")]
chk("homicidios entre matrices", f"{es(100 * min(cai), 0)} % a {es(100 * max(cai), 0)} %", "34 % a 41 %", 0, "La caída relativa en los homicidios va de 34 % a 41 % según la matriz")
NOM_A2 = {"todas": "Todas las causas", "604": "Alzheimer, Parkinson y otras", "201": "Cáncer de estómago", "613": "Malformaciones congénitas", "204": "Cáncer de pulmón",
          "406": "Otras afecciones perinatales", "108": "Neumonía e IRA", "210": "Cáncer de próstata", "212": "Leucemia", "302": "Enf. hipertensivas",
          "505": "Ahogamiento", "404": "Dificultad respiratoria del RN", "103": "Enf. por vectores y rabia", "512": "Homicidios"}
for c, nom in NOM_A2.items():
    r = fila(T4, causa=c)
    linea = f"| {nom} | " + " | ".join(f"{es(r['I_y_' + w])} → {es(r['I_SGQ_' + w])}" for w in ("reina", "knn5", "knn8", "invdist")) + " |"
    chk(f"tabla A2 {c}", linea, linea, 0, linea)
chk("sin pandemia", f"{es(C['todas']['moran_y']['I'])} a {es(C['todas']['sin_pandemia']['moran_y'])} y de {es(C['todas']['I_resid']['SGQ']['reina'])} a {es(C['todas']['sin_pandemia']['moran_res'])}",
    "0,29 a 0,30 y de 0,18 a 0,20", 0, "pasan de 0,29 a 0,30 y de 0,18 a 0,20")
chk("Conley de MCO", (T6["se_conley100"] / T6["se_ols"]).median(), "1,33", 2, "1,33 veces los clásicos")
chk("GMM", max(abs(C[c]["gmm_sdem"]["lam"] - C[c]["modelos"]["SDEM"]["lam"]) for c in C), "0,02", 2, "difiere en 0,02 o menos del de máxima verosimilitud")
ET = json.loads((RAIZ / "data" / "raw" / "externos" / "63" / "etnia2018_mpio.json").read_text())["features"]
muns = set(COV["mun"].to_list())
nul = lambda k: sum(1 for x in ET if int(x["attributes"]["U_MPIO"]) in muns and x["attributes"][k] is None)
chk("celdas vacías", f"{nul('PARTICIPACION_PORCENTUAL_PA1_GRP_ETNIC_1')} y {nul('PARTICIPACION_PORCENTUAL_PA1_GRP_ETNIC_5')}", "152 y 30", 0, "se leen como cero (152 y 30 municipios)")

if fallas:
    print(f"{len(fallas)} fallas de {len(fallas) + n_ok} comprobaciones:")
    for f in fallas:
        print("  ·", f)
    sys.exit(1)
print(f"{n_ok} comprobaciones.")
print("Todo coincide.")
