"""Recálculo independiente de las cifras del resumen del manuscrito.

No usa el paquete `agrupamiento`: lee data/intermedios/ con pandas y recalcula con esda, libpysal, statsmodels y código
propio (parte del recálculo del árbitro). Cifras del resumen que se comprueban:

  · 2,23 millones de defunciones, 2018–2025, 1.120 municipios;
  · la I de Moran inicial de las series (0,06 a 0,62) y las 12 series con I > 0,1;
  · covariables placebo con la misma estructura espacial reducen la I entre 10 % y 22 % (en promedio, por serie);
  · las reales la reducen entre 13 y 33 puntos más en 11 de 12 series, pero no en las enfermedades hipertensivas;
  · la calidad del registro aporta más que el placebo en seis series;
  · sobreviven 14 de 224 asociaciones (SLX ponderado con errores agrupados recalculado con statsmodels; el lado BYM2 se
    toma de output/bym2.json, cuyo estimador se prueba contra casos conocidos en test_estimadores.py).

El placebo se recalcula con otra semilla, otra calibración (bisección) y 100 réplicas: se admite una tolerancia de
1,5 puntos en los promedios y se exige que las conclusiones (qué series superan al placebo) sean las mismas.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

RAIZ = Path(__file__).resolve().parents[1]
INT = RAIZ / "data" / "intermedios"
OUT = RAIZ / "output"
COVS = ["nbi", "rural", "ldens", "indigena", "afro", "internet", "superior", "altitud", "ldist_capital", "ldist_uci",
        "lcamas", "pdet", "sin_asistencia", "basura", "casa", "completitud"]
BL = {"S": COVS[:7], "G": COVS[7:12], "Q": COVS[12:]}


@pytest.fixture(scope="module")
def base():
    import libpysal
    cov = pd.read_parquet(INT / "covariables.parquet").sort_values("mun").reset_index(drop=True)
    oe = pd.read_parquet(INT / "observados_esperados.parquet")
    adj = pd.read_csv(INT / "pesos_reina.csv")
    ids = cov["mun"].tolist()
    nb = {m: [] for m in ids}
    for i, j in adj.itertuples(index=False):
        nb[i].append(j)
    w = libpysal.weights.W(nb, ids=ids, silence_warnings=True); w.transform = "r"
    Wd = w.full()[0]
    Y = {}
    for c, g in oe.groupby("causa"):
        g = g.set_index("mun").loc[ids]
        O, E = g["O"].to_numpy(float), g["E"].to_numpy()
        r = O / E; m = O.sum() / E.sum()
        s2 = max((E * (r - m) ** 2).sum() / E.sum() - m / E.mean(), 0)
        Y[c] = np.log(m + s2 / (s2 + m / E) * (r - m))
    Z = {v: (cov[v] - cov[v].mean()) / cov[v].std(ddof=0) if v != "pdet" else cov[v].astype(float) for v in COVS}
    X = pd.DataFrame(Z).to_numpy()
    return {"cov": cov, "oe": oe, "w": w, "Wd": Wd, "Y": Y, "X": X}


def _I(e, Wd):
    e = e - e.mean(); return len(e) / Wd.sum() * (e @ Wd @ e) / (e @ e)


def _res(y, X):
    Zm = np.column_stack([np.ones(len(y)), X]) if X.size else np.ones((len(y), 1))
    return y - Zm @ np.linalg.lstsq(Zm, y, rcond=None)[0]


def test_tamano_muestra(base):
    oe = base["oe"]
    assert round(oe.loc[oe.causa == "todas", "O"].sum() / 1e6, 2) == 2.23
    assert base["cov"].shape[0] == 1120


def test_moran_inicial(base):
    import esda
    I = {c: esda.Moran(y, base["w"], permutations=0).I for c, y in base["Y"].items()}
    assert round(min(I.values()), 2) == 0.06 and round(max(I.values()), 2) == 0.62
    assert sum(v > 0.1 for v in I.values()) == 12


def _calibrar(I_obj, Wd, rng, n):
    lo, hi = 0.0, 0.98
    for _ in range(18):
        r = (lo + hi) / 2
        A = np.linalg.inv(np.eye(n) - r * Wd)
        Im = np.mean([_I(A @ rng.standard_normal(n), Wd) for _ in range(20)])
        lo, hi = (r, hi) if Im < I_obj else (lo, r)
    return (lo + hi) / 2


@pytest.fixture(scope="module")
def placebo(base):
    rng = np.random.default_rng(99)
    Wd, X, cov = base["Wd"], base["X"], base["cov"]; n = Wd.shape[0]
    Iv = [_I(X[:, j], Wd) for j in range(X.shape[1])]
    A = [np.linalg.inv(np.eye(n) - _calibrar(Iv[j], Wd, rng, n) * Wd) for j in range(X.shape[1])]
    doce = [c for c, y in base["Y"].items() if _I(y, Wd) > 0.1]
    pos = {b: [COVS.index(v) for v in vs] for b, vs in BL.items()}
    sub = lambda s: np.column_stack([X[:, pos[b]] for b in s]) if s else np.empty((n, 0))
    out = {}
    nrep = 100
    sims = []
    for _ in range(nrep):
        N = np.column_stack([A[j] @ rng.standard_normal(n) for j in range(X.shape[1])])
        sims.append((N - N.mean(0)) / N.std(0))
    for c in doce:
        y = base["Y"][c]; I0 = _I(y, Wd)
        Ireal = {"".join(s) or "0": _I(_res(y, sub(s)), Wd) for r_ in range(4) for s in itertools.combinations("SGQ", r_)}
        caida_real = 1 - Ireal["SGQ"] / I0
        caidas = np.array([1 - _I(_res(y, N), Wd) / I0 for N in sims])
        # exceso de Shapley del bloque Q sobre un bloque placebo del mismo tamaño
        phi_real = sum(math.factorial(len(s)) * math.factorial(2 - len(s)) / 6 *
                       (Ireal["".join(x for x in "SG" if x in s) or "0"] - Ireal["".join(x for x in "SGQ" if x in s or x == "Q")])
                       for r_ in range(3) for s in itertools.combinations("SG", r_))
        phi_pl = []
        for N in sims:
            Nq = N[:, pos["Q"]]
            phi_pl.append(sum(math.factorial(len(s)) * math.factorial(2 - len(s)) / 6 *
                              (Ireal["".join(x for x in "SG" if x in s) or "0"] - _I(_res(y, np.column_stack([sub(s), Nq]) if s else Nq), Wd))
                              for r_ in range(3) for s in itertools.combinations("SG", r_)))
        phi_pl = np.array(phi_pl)
        out[c] = {"real": caida_real, "placebo": caidas.mean(), "p": (np.sum(caidas >= caida_real) + 1) / (nrep + 1),
                  "q_lo": (phi_real - np.quantile(phi_pl, 0.975)) / I0}
    return out


def test_placebo_entre_10_y_22(placebo):
    m = [v["placebo"] for v in placebo.values()]
    assert abs(100 * min(m) - 10.2) < 1.5 and abs(100 * max(m) - 21.8) < 1.5


def test_exceso_13_a_33_en_11_de_12(placebo):
    sup = {c: v for c, v in placebo.items() if v["p"] <= 0.05}
    assert len(sup) == 11 and "302" not in sup
    ex = [100 * (v["real"] - v["placebo"]) for v in sup.values()]
    assert abs(min(ex) - 12.5) < 1.5 and abs(max(ex) - 32.8) < 1.5


def test_calidad_del_registro_supera_al_placebo_en_seis(placebo):
    q = sorted(c for c, v in placebo.items() if v["q_lo"] > 0)
    assert q == ["103", "201", "212", "404", "512", "todas"]


def test_sobreviven_14_de_224(base):
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests
    from scipy.stats import norm
    cov, Wd, X = base["cov"], base["Wd"], base["X"]
    Xs = X.copy(); Xs[:, COVS.index("pdet")] = cov["pdet"].to_numpy(float)
    Zm = sm.add_constant(np.column_stack([Xs, Wd @ Xs]))
    k = len(COVS); causas = sorted(base["Y"])
    filas = []
    by = json.loads((OUT / "bym2.json").read_text())["causas"]
    oe = base["oe"]
    for c in causas:
        E = oe[oe.causa == c].set_index("mun").loc[cov["mun"], "E"].to_numpy()
        m = sm.WLS(base["Y"][c], Zm, weights=E / E.mean()).fit(cov_type="cluster", cov_kwds={"groups": cov["dp"].to_numpy()})
        V = m.cov_params(); b = m.params
        bb = np.array(by[c]["completo"]["beta"]); Vb = np.array(by[c]["completo"]["cov_beta"])
        for j in range(k):
            t = b[1 + j] + b[1 + k + j]; s = math.sqrt(V[1 + j, 1 + j] + V[1 + k + j, 1 + k + j] + 2 * V[1 + j, 1 + k + j])
            tb = bb[1 + j] + bb[1 + k + j]; sb = math.sqrt(Vb[1 + j, 1 + j] + Vb[1 + k + j, 1 + k + j] + 2 * Vb[1 + j, 1 + k + j])
            filas.append((c, COVS[j], t, 2 * (1 - norm.cdf(abs(t / s))), tb, 2 * (1 - norm.cdf(abs(tb / sb)))))
    T = pd.DataFrame(filas, columns=["c", "v", "t", "p", "tb", "pb"])
    T["q"] = multipletests(T["p"], method="fdr_bh")[1]; T["qb"] = multipletests(T["pb"], method="fdr_bh")[1]
    res = json.loads((OUT / "resultados.json").read_text())["causas"]
    T["sdem"] = [res[c]["impactos"]["SDEM"]["total"]["media"][COVS.index(v)] for c, v in zip(T.c, T.v)]
    sob = (T.q < 0.05) & (T.qb < 0.05) & (np.sign(T.t) == np.sign(T.tb)) & (np.sign(T.sdem) == np.sign(T.t))
    assert len(T) == 224 and int((T.q < 0.05).sum()) == 52 and int(sob.sum()) == 14
