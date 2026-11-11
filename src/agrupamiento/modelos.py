"""Etapas de estimación del nivel A (todas parten de data/intermedios/).

  etapa_modelos()        Moran, MCO por bloques, Shapley, LM, SAR/SEM/SDM/SDEM/SLX, efectos, Conley, GMM, otras
                         matrices, regímenes, años sin pandemia                     -> output/resultados.json
  etapa_efectos_fijos()  SDM y SDEM con efectos fijos de departamento              -> resultados.json["causas"][c]["fe_espacial"]
  etapa_revision()       placebo calibrado, numerador y denominador, E[I], SLX ponderado y sin pesos con errores
                         agrupados y de Conley, BH, códigos basura sin I10–I15     -> output/revision.json
  etapa_bym2()           Poisson-BYM2 con rezagos de las covariables               -> output/bym2.json
  etapa_panel()          panel anual con rezago espacio-temporal                   -> output/panel.json
  etapa_departamental()  I de Moran departamental de 2022 de las 14 series          -> output/departamental.json
"""
from __future__ import annotations

import itertools
import json
import math
import time

import numpy as np
import polars as pl

from . import pesos, rutas
from .constantes import (BLOQUES, CAUSAS, CLAVES, COVS, SEMILLA, SEMILLA_PLACEBO, clave)
from .estimadores import (bayes_empirico, bh, bym2_laplace, conley_se, escala_icar, esperanza_moran_h0, feols_dos_vias,
                          impactos_sdem, impactos_sdm, moran, moran_columnas, numerador_denominador, ols_resid, shapley,
                          slx_inferencia, zstd)


def _json(obj, ruta):
    ruta.parent.mkdir(parents=True, exist_ok=True)
    ruta.write_text(json.dumps(obj, ensure_ascii=False, default=float))


def cargar():
    cov = pl.read_parquet(rutas.INTER / "covariables.parquet").sort("mun")
    oe = pl.read_parquet(rutas.INTER / "observados_esperados.parquet")
    return cov, oe


def y_causa(cov, oe, c):
    d = cov.select("mun").join(oe.filter(pl.col("causa") == c), on="mun", how="left", maintain_order="left")
    O, E = d["O"].to_numpy().astype(float), d["E"].to_numpy()
    th, eb = bayes_empirico(O, E)
    return np.log(th), O, E, eb


def disenio(cov: pl.DataFrame, vars_, fe=False):
    if fe:   # la completitud es departamental: la absorben los efectos fijos
        vars_ = [v for v in vars_ if v != "completitud"]
    X = np.column_stack([zstd(cov[v].to_numpy()) if v != "pdet" else cov[v].to_numpy().astype(float) for v in vars_]) if vars_ else np.empty((cov.height, 0))
    nombres = list(vars_)
    if fe:
        dps = sorted(set(cov["dp"].to_list()))
        D = np.column_stack([(cov["dp"].to_numpy() == d).astype(float) for d in dps[1:]])
        X = np.column_stack([X, D]) if X.size else D
        nombres += [f"dp_{d}" for d in dps[1:]]
    return X, nombres


def ajustar_causa(k, y, O, E, cov, W, Wd, eig, rng):
    import esda
    import spreg
    from scipy.stats import chi2
    res = {"causa": clave(k), "etiqueta": CAUSAS[k][0], "O": int(O.sum()), "E": float(E.sum())}
    wq = W["reina"]
    np.random.seed(SEMILLA)                                   # permutaciones de esda: deterministas
    mi = esda.Moran(y, wq, permutations=999)
    mr = esda.Moran_Rate(O.astype(float), E, wq, adjusted=True, permutations=999)
    res["moran_y"] = {"I": float(mi.I), "p_sim": float(mi.p_sim), "z": float(mi.z_norm)}
    res["moran_ar"] = {"I": float(mr.I), "p_sim": float(mr.p_sim)}
    Ires = {}
    for r_ in range(4):
        for sub in itertools.combinations(list(BLOQUES), r_):
            X, _ = disenio(cov, [v for b_ in sub for v in BLOQUES[b_]])
            e, _ = ols_resid(y, X)
            Ires["".join(sub) or "0"] = {w: moran(e, Wd[w]) for w in Wd}
    Xf, nf = disenio(cov, COVS, fe=True)
    Ires["SGQ+FE"] = {w: moran(ols_resid(y, Xf)[0], Wd[w]) for w in Wd}
    Ires["FE"] = {w: moran(ols_resid(y, disenio(cov, [], fe=True)[0])[0], Wd[w]) for w in Wd}
    res["I_resid"] = Ires
    res["shapley"] = {w: shapley({s: v[w] for s, v in Ires.items()}) for w in Wd}
    X, nombres = disenio(cov, COVS)
    ols = spreg.OLS(y[:, None], X, w=wq, spat_diag=True, moran=True, name_x=nombres)
    res["ols"] = {"beta": np.asarray(ols.betas).ravel().tolist(), "se": np.asarray(ols.std_err).ravel().tolist(),
                  "r2": float(ols.r2), "ar2": float(ols.ar2), "aic": float(ols.aic), "ll": float(ols.logll),
                  "moran_res": [float(x) for x in ols.moran_res],
                  "lm_error": [float(x) for x in ols.lm_error], "lm_lag": [float(x) for x in ols.lm_lag],
                  "rlm_error": [float(x) for x in ols.rlm_error], "rlm_lag": [float(x) for x in ols.rlm_lag],
                  "lm_sarma": [float(x) for x in ols.lm_sarma], "nombres": nombres,
                  "vif": np.diag(np.linalg.inv(np.corrcoef(X.T))).tolist()}
    ols_fe = spreg.OLS(y[:, None], Xf, w=wq, spat_diag=True, moran=True, nonspat_diag=False)
    res["ols_fe"] = {"moran_res": [float(x) for x in ols_fe.moran_res], "lm_error": [float(x) for x in ols_fe.lm_error],
                     "lm_lag": [float(x) for x in ols_fe.lm_lag], "rlm_error": [float(x) for x in ols_fe.rlm_error],
                     "rlm_lag": [float(x) for x in ols_fe.rlm_lag], "ar2": float(ols_fe.ar2), "nombres": nf[:len(COVS) - 1],
                     "beta": np.asarray(ols_fe.betas).ravel()[:len(COVS)].tolist(), "se": np.asarray(ols_fe.std_err).ravel()[:len(COVS)].tolist()}
    _, ses = conley_se(y, X, cov["lat"].to_numpy(), cov["lon"].to_numpy())
    res["conley"] = {str(c): v.tolist() for c, v in ses.items()}
    nk = len(COVS); modelos = {}

    def guardar(nom, m, tipo):
        d = {"beta": np.asarray(m.betas).ravel().tolist(), "se": np.asarray(m.std_err).ravel().tolist(),
             "ll": float(m.logll), "aic": float(m.aic), "nombres": list(m.name_x) + (["rho"] if tipo == "lag" else ["lambda"])}
        d["rho" if tipo == "lag" else "lam"] = float(m.rho if tipo == "lag" else m.lam)
        modelos[nom] = d
        return m
    sar = guardar("SAR", spreg.ML_Lag(y[:, None], X, wq, method="ord", spat_impacts=None, spat_diag=False, name_x=nombres), "lag")
    sem = guardar("SEM", spreg.ML_Error(y[:, None], X, wq, method="ord", name_x=nombres), "err")
    sdm = guardar("SDM", spreg.ML_Lag(y[:, None], X, wq, slx_lags=1, method="ord", spat_impacts=None, spat_diag=False, name_x=nombres), "lag")
    sdem = guardar("SDEM", spreg.ML_Error(y[:, None], X, wq, slx_lags=1, method="ord", name_x=nombres), "err")
    slx = spreg.OLS(y[:, None], X, w=wq, slx_lags=1, name_x=nombres)
    modelos["SLX"] = {"beta": np.asarray(slx.betas).ravel().tolist(), "se": np.asarray(slx.std_err).ravel().tolist(),
                      "ll": float(slx.logll), "aic": float(slx.aic)}
    modelos["OLS"] = {"ll": float(ols.logll), "aic": float(ols.aic)}

    def lr(a, b_, df):
        s = 2 * (modelos[a]["ll"] - modelos[b_]["ll"]); return [float(s), float(chi2.sf(s, df))]
    res["lr"] = {"SDM_vs_SAR": lr("SDM", "SAR", nk), "SDM_vs_SEM": lr("SDM", "SEM", nk), "SDM_vs_SLX": lr("SDM", "SLX", 1),
                 "SDEM_vs_SEM": lr("SDEM", "SEM", nk), "SDEM_vs_SLX": lr("SDEM", "SLX", 1), "SLX_vs_OLS": lr("SLX", "OLS", nk),
                 "SEM_vs_OLS": lr("SEM", "OLS", 1), "SAR_vs_OLS": lr("SAR", "OLS", 1)}
    res["modelos"] = modelos
    res["moran_res_sdem"] = moran(np.asarray(sdem.e_filtered).ravel(), Wd["reina"])
    res["moran_res_sem"] = moran(np.asarray(sem.e_filtered).ravel(), Wd["reina"])
    res["moran_res_sdm"] = moran(np.asarray(sdm.u).ravel(), Wd["reina"])
    imp = {"SDM": impactos_sdm(sdm, nk, eig["reina"], 2000, rng), "SAR": impactos_sdm(sar, nk, eig["reina"], 2000, rng),
           "SDEM": impactos_sdem(sdem, nk)}
    res["impactos"] = {m: {t: {kk: np.asarray(vv).tolist() for kk, vv in d.items()} for t, d in im.items()} for m, im in imp.items()}
    gh = spreg.GM_Error_Het(y[:, None], X, wq, slx_lags=1, name_x=nombres)
    b_ = np.asarray(gh.betas).ravel()
    res["gmm_sdem"] = {"beta": b_.tolist(), "se": np.asarray(gh.std_err).ravel().tolist(), "lam": float(b_[-1])}
    otras = {}
    for wn in ("knn5", "knn8", "invdist"):
        w = W[wn]
        o_ = spreg.OLS(y[:, None], X, w=w, spat_diag=True, moran=True)
        sd = spreg.ML_Error(y[:, None], X, w, slx_lags=1, method="ord")
        sm = spreg.ML_Lag(y[:, None], X, w, slx_lags=1, method="ord", spat_impacts=None, spat_diag=False)
        otras[wn] = {"moran_y": float(esda.Moran(y, w, permutations=0).I), "moran_res": [float(x) for x in o_.moran_res],
                     "rlm_error": [float(x) for x in o_.rlm_error], "rlm_lag": [float(x) for x in o_.rlm_lag],
                     "aic_sdem": float(sd.aic), "aic_sdm": float(sm.aic), "lam_sdem": float(sd.lam), "rho_sdm": float(sm.rho),
                     "imp_sdem": {t: {kk: np.asarray(vv).tolist() for kk, vv in d.items()} for t, d in impactos_sdem(sd, nk).items()},
                     "imp_sdm_total": {kk: np.asarray(vv).tolist() for kk, vv in impactos_sdm(sm, nk, eig[wn], 1000, rng)["total"].items()}}
    res["otras_W"] = otras
    rg = spreg.OLS_Regimes(y[:, None], X, regimes=cov["macro"].to_list(), constant_regi="many", cols2regi="all",
                           regime_err_sep=False, name_x=nombres)
    res["regimenes"] = {"chow_global": [float(rg.chow.joint[0]), float(rg.chow.joint[1])],
                        "chow_vars": {nombres[i]: [float(rg.chow.regi[i + 1, 0]), float(rg.chow.regi[i + 1, 1])] for i in range(nk)}}
    return res


def etapa_modelos():
    rng = np.random.default_rng(SEMILLA)
    cov, oe = cargar()
    oe_sp = pl.read_parquet(rutas.INTER / "observados_esperados_sin_pandemia.parquet")
    W, winfo = pesos.todas(cov)
    Wd = {k: w.full()[0] for k, w in W.items()}
    eig = {k: np.linalg.eigvals(v) for k, v in Wd.items()}
    salida = {"meta": {"n": cov.height, "W": winfo, "semilla": SEMILLA, "covariables": COVS, "bloques": BLOQUES}, "causas": {}}
    for k in CLAVES:
        t0 = time.time()
        y, O, E, eb = y_causa(cov, oe, clave(k))
        r = ajustar_causa(k, y, O, E, cov, W, Wd, eig, rng)
        r["eb"] = {kk: v for kk, v in eb.items() if kk != "C"}
        r["ceros"] = int((O == 0).sum())
        y2, *_ = y_causa(cov, oe_sp, clave(k))
        X, _ = disenio(cov, COVS)
        e2, b2 = ols_resid(y2, X)
        r["sin_pandemia"] = {"O": int(oe_sp.filter(pl.col("causa") == clave(k))["O"].sum()),
                             "moran_y": moran(y2, Wd["reina"]), "moran_res": moran(e2, Wd["reina"]), "beta": b2.tolist()}
        salida["causas"][clave(k)] = r
        print(f"  {clave(k):>6} I={r['moran_y']['I']:.3f} Ires={r['I_resid']['SGQ']['reina']:.3f} {time.time() - t0:.0f}s", flush=True)
    _json(salida, rutas.OUT / "resultados.json")


def etapa_efectos_fijos():
    import spreg
    cov, oe = cargar()
    W, _ = pesos.todas(cov)
    wq = W["reina"]
    R = json.loads((rutas.OUT / "resultados.json").read_text())
    for k in CLAVES:
        y, *_ = y_causa(cov, oe, clave(k))
        Xf, _ = disenio(cov, COVS, fe=True)
        kc = len([v for v in COVS if v != "completitud"])
        sv = [True] * kc + [False] * (Xf.shape[1] - kc)
        sdm = spreg.ML_Lag(y[:, None], Xf, wq, slx_lags=1, slx_vars=sv, method="ord", spat_impacts=None, spat_diag=False)
        sdem = spreg.ML_Error(y[:, None], Xf, wq, slx_lags=1, slx_vars=sv, method="ord")
        sem = spreg.ML_Error(y[:, None], Xf, wq, method="ord")
        sar = spreg.ML_Lag(y[:, None], Xf, wq, method="ord", spat_impacts=None, spat_diag=False)
        ols = spreg.OLS(y[:, None], Xf, w=wq, slx_lags=1, slx_vars=sv, spat_diag=True, moran=True, nonspat_diag=False)
        se = lambda m: float(np.asarray(m.std_err).ravel()[-1])
        R["causas"][clave(k)]["fe_espacial"] = {
            "rho_sdm": float(sdm.rho), "se_rho_sdm": se(sdm), "aic_sdm": float(sdm.aic), "ll_sdm": float(sdm.logll),
            "lam_sdem": float(sdem.lam), "se_lam_sdem": se(sdem), "aic_sdem": float(sdem.aic), "ll_sdem": float(sdem.logll),
            "lam_sem": float(sem.lam), "rho_sar": float(sar.rho), "aic_sem": float(sem.aic), "aic_sar": float(sar.aic),
            "slx_moran_res": [float(x) for x in ols.moran_res], "slx_rlm_error": [float(x) for x in ols.rlm_error],
            "slx_rlm_lag": [float(x) for x in ols.rlm_lag]}
        print(f"  EF {clave(k):>6} rho={sdm.rho:.3f} lam={sdem.lam:.3f}", flush=True)
    _json(R, rutas.OUT / "resultados.json")


def etapa_revision(nrep: int = 200):
    """Línea base placebo (ecuación 5), numerador y denominador de la I, E[I] bajo H0, SLX con y sin pesos de precisión
    con errores agrupados y de Conley, BH conjunto sobre el efecto total, códigos basura sin I10–I15."""
    from scipy.stats import norm
    rng = np.random.default_rng(SEMILLA_PLACEBO)
    cov, oe = cargar()
    R = json.loads((rutas.OUT / "resultados.json").read_text())["causas"]
    W, _ = pesos.todas(cov)
    Wd = W["reina"].full()[0]; n = Wd.shape[0]
    Ys, Es = {}, {}
    for k in CLAVES:
        y, O, E, _ = y_causa(cov, oe, clave(k))
        Ys[clave(k)] = y; Es[clave(k)] = E
    cl = list(Ys); Y = np.column_stack([Ys[c] for c in cl])

    def resid(Xm):
        Z = np.column_stack([np.ones(n), Xm]) if Xm.size else np.ones((n, 1))
        return Y - Z @ np.linalg.lstsq(Z, Y, rcond=None)[0]

    Xb = {b: disenio(cov, BLOQUES[b])[0] for b in BLOQUES}
    sub = lambda s: np.column_stack([Xb[b] for b in s]) if s else np.empty((n, 0))
    out = {"meta": {"replicas": nrep, "semilla": SEMILLA_PLACEBO}}
    out["EI_H0"] = {nom: esperanza_moran_h0(Xm, Wd) for nom, Xm in (("0", np.empty((n, 0))), ("SGQ", sub("SGQ")),
                                                                     ("FE", disenio(cov, [], fe=True)[0]),
                                                                     ("SGQ+FE", disenio(cov, COVS, fe=True)[0]))}
    reales = {"".join(s) or "0": resid(sub(s)) for r_ in range(4) for s in itertools.combinations("SGQ", r_)}
    Ireal = {s: moran_columnas(E_, Wd) for s, E_ in reales.items()}
    nd0 = numerador_denominador(reales["0"], Wd); nd1 = numerador_denominador(reales["SGQ"], Wd)
    efe = resid(disenio(cov, COVS, fe=True)[0]); ndfe = numerador_denominador(efe, Wd)
    # placebo calibrado: ruido autorregresivo con la I de Moran de cada covariable real
    Iv = {v: float(moran_columnas(zstd(cov[v].to_numpy().astype(float))[:, None], Wd)[0]) for v in COVS}
    grid = np.linspace(0, 0.97, 25)
    Ig = []
    for rho in grid:
        Ig.append(float(np.mean(moran_columnas(np.linalg.inv(np.eye(n) - rho * Wd) @ rng.standard_normal((n, 40)), Wd))))
    Ig = np.maximum.accumulate(np.array(Ig))
    rho_v = {v: float(np.interp(Iv[v], Ig, grid)) for v in COVS}
    Ainvs = {r: np.linalg.inv(np.eye(n) - r * Wd) for r in sorted(set(round(r, 3) for r in rho_v.values()))}
    pos = {b: [COVS.index(v) for v in BLOQUES[b]] for b in BLOQUES}
    caida_tot = np.zeros((nrep, len(cl))); phi_pl = {b: np.zeros((nrep, len(cl))) for b in BLOQUES}
    caida_blk = {b: np.zeros((nrep, len(cl))) for b in BLOQUES}
    I0 = Ireal["0"]
    for rep in range(nrep):
        N = np.column_stack([zstd(Ainvs[round(rho_v[v], 3)] @ rng.standard_normal(n)) for v in COVS])
        Nb = {b: N[:, pos[b]] for b in BLOQUES}
        caida_tot[rep] = 1 - moran_columnas(resid(N), Wd) / I0
        for b in BLOQUES:
            caida_blk[b][rep] = 1 - moran_columnas(resid(Nb[b]), Wd) / I0
            otros = [x for x in "SGQ" if x != b]; acc = np.zeros(len(cl))
            for r_ in range(3):
                for s in itertools.combinations(otros, r_):
                    peso = math.factorial(len(s)) * math.factorial(2 - len(s)) / 6
                    s_sin = "".join(x for x in "SGQ" if x in s) or "0"
                    Ic = moran_columnas(resid(np.column_stack([sub(s), Nb[b]]) if s else Nb[b]), Wd)
                    acc += peso * (Ireal[s_sin] - Ic)
            phi_pl[b][rep] = acc
    res = {}
    for j, c in enumerate(cl):
        real_tot = 1 - Ireal["SGQ"][j] / I0[j]
        d = {"I0": float(I0[j]), "I_SGQ": float(Ireal["SGQ"][j]), "I_SGQFE": float(moran_columnas(efe[:, [j]], Wd)[0]),
             "caida_real": float(real_tot), "caida_placebo_media": float(caida_tot[:, j].mean()),
             "caida_placebo_p95": float(np.quantile(caida_tot[:, j], 0.95)),
             "caida_exceso": float(real_tot - caida_tot[:, j].mean()),
             "p_placebo": float((np.sum(caida_tot[:, j] >= real_tot) + 1) / (nrep + 1)),
             "I_SGQ_placebo_p05_p50_p95": np.quantile(I0[j] * (1 - caida_tot[:, j]), [0.05, 0.5, 0.95]).tolist(),
             "num0": float(nd0[0][j]), "den0": float(nd0[1][j]), "num_SGQ": float(nd1[0][j]), "den_SGQ": float(nd1[1][j]),
             "num_FE": float(ndfe[0][j]), "den_FE": float(ndfe[1][j]),
             "caida_num": float(1 - nd1[0][j] / nd0[0][j]), "caida_den": float(1 - nd1[1][j] / nd0[1][j]), "bloques": {}}
        for b in BLOQUES:
            phr = R[c]["shapley"]["reina"][b]; pl_ = phi_pl[b][:, j]
            d["bloques"][b] = {"n_vars": len(BLOQUES[b]), "shapley_real": float(phr / I0[j]),
                               "shapley_placebo_media": float(pl_.mean() / I0[j]),
                               "shapley_exceso": float((phr - pl_.mean()) / I0[j]),
                               "shapley_exceso_ic": [float((phr - np.quantile(pl_, 0.975)) / I0[j]), float((phr - np.quantile(pl_, 0.025)) / I0[j])],
                               "p_placebo": float((np.sum(pl_ >= phr) + 1) / (nrep + 1)),
                               "caida_sola_real": float(1 - Ireal[b][j] / I0[j]),
                               "caida_sola_placebo_media": float(caida_blk[b][:, j].mean()),
                               "caida_sola_placebo_p95": float(np.quantile(caida_blk[b][:, j], 0.95))}
        Ir = {s: Ireal[s][j] for s in Ireal}
        d["no_monotono"] = [[a, bb] for a, bb in (("0", "S"), ("0", "G"), ("0", "Q"), ("S", "SG"), ("S", "SQ"), ("G", "SG"),
                                                   ("G", "GQ"), ("Q", "SQ"), ("Q", "GQ"), ("SG", "SGQ"), ("SQ", "SGQ"), ("GQ", "SGQ"))
                            if Ir[bb] > Ir[a] + 1e-9]
        res[c] = d
    out["causas"] = res
    out["I_covariables"] = Iv; out["rho_placebo"] = rho_v
    X, _ = disenio(cov, COVS)
    X2, _ = disenio(cov, [v if v != "basura" else "basura_sin_hta" for v in COVS])
    g = cov["dp"].to_numpy(); la, lo = cov["lat"].to_numpy(), cov["lon"].to_numpy()
    slx = {c: {"MCO": slx_inferencia(Ys[c], X, Wd, np.ones(n), g, la, lo),
               "MCP": slx_inferencia(Ys[c], X, Wd, Es[c] / Es[c].mean(), g, la, lo),
               "MCP_sin_hta": slx_inferencia(Ys[c], X2, Wd, Es[c] / Es[c].mean(), g, la, lo)} for c in cl}
    for mod in ("MCO", "MCP"):
        for se_ in ("cl", "conley"):
            p = []
            for c in cl:
                t = np.array(slx[c][mod]["tot"]); s = np.array(slx[c][mod][se_]["se_tot"])
                p += list(2 * (1 - norm.cdf(np.abs(t / s))))
            q = bh(p)
            for i_, c in enumerate(cl):
                slx[c][mod][se_]["q_tot"] = q[i_ * len(COVS):(i_ + 1) * len(COVS)].tolist()
                slx[c][mod][se_]["p_tot"] = p[i_ * len(COVS):(i_ + 1) * len(COVS)]
    out["slx"] = slx
    _json(out, rutas.OUT / "revision.json")


def etapa_bym2(solo=None):
    """BYM2 escalado por Laplace: nulo y con S, G, Q y sus rezagos espaciales; para las hipertensivas, además, con los
    códigos basura sin I10–I15. Guarda después de cada causa (se puede interrumpir y seguir)."""
    cov, oe = cargar()
    muns = cov["mun"].to_list()
    W, _ = pesos.todas(cov)
    Wd = W["reina"].full()[0]
    ev, U, c_esc = escala_icar(pesos.simetrica_binaria(muns))
    X, _ = disenio(cov, COVS); XW = np.column_stack([X, Wd @ X])
    ruta = rutas.OUT / "bym2.json"
    salida = json.loads(ruta.read_text()) if (ruta.exists() and solo != "todo") else {"meta": {"c_escala": c_esc}, "causas": {}}
    for k in CLAVES:
        c = clave(k)
        if c in salida["causas"]:
            continue
        t0 = time.time()
        d = cov.select("mun").join(oe.filter(pl.col("causa") == c), on="mun", how="left", maintain_order="left")
        O, E = d["O"].to_numpy().astype(float), d["E"].to_numpy()
        m0 = bym2_laplace(O, E, np.empty((len(O), 0)), ev, U, c_esc)
        m1 = bym2_laplace(O, E, XW, ev, U, c_esc, theta0=m0["theta"])
        r = {"nulo": {kk: v for kk, v in m0.items() if kk != "u"}, "completo": {kk: v for kk, v in m1.items() if kk != "u"},
             "moran_u_nulo": moran(np.array(m0["u"]), Wd), "moran_u_completo": moran(np.array(m1["u"]), Wd)}
        if c == "302":
            X2, _ = disenio(cov, [v if v != "basura" else "basura_sin_hta" for v in COVS])
            m2 = bym2_laplace(O, E, np.column_stack([X2, Wd @ X2]), ev, U, c_esc, theta0=m1["theta"])
            r["sin_hta"] = {kk: v for kk, v in m2.items() if kk != "u"}
        r["seg"] = time.time() - t0
        salida["causas"][c] = r
        _json(salida, ruta)
        print(f"  BYM2 {c:>6} σ² {m0['sigma2']:.3f}→{m1['sigma2']:.3f} φ {m0['phi']:.3f}→{m1['phi']:.3f} {r['seg']:.0f}s", flush=True)


def etapa_panel():
    """Panel anual 2018–2025 de log SMR-EB con rezago propio y de los vecinos del año anterior (ecuación 8)."""
    cov, _ = cargar()
    W, _ = pesos.todas(cov)
    Wd = W["reina"].full()[0]
    muns = cov["mun"].to_list()
    anual = pl.read_parquet(rutas.INTER / "observados_esperados_anual.parquet")
    out = {}
    for k in CLAVES:
        c = clave(k); filas = []
        for a in sorted(anual["anio"].unique().to_list()):
            d = pl.DataFrame({"mun": muns}).join(anual.filter((pl.col("causa") == c) & (pl.col("anio") == a)), on="mun", how="left", maintain_order="left")
            O_, E_ = d["O"].fill_null(0).to_numpy().astype(float), d["E"].fill_null(0).to_numpy()
            ok = E_ > 0                       # Mapiripana (94663) no tiene población proyectada desde 2020
            y = np.zeros(len(O_)); y[ok] = np.log(bayes_empirico(O_[ok], E_[ok])[0])
            filas.append(pl.DataFrame({"mun": muns, "anio": a, "y": np.where(ok, y, np.nan), "wy": Wd @ y}))
        P = pl.concat(filas).sort("mun", "anio").with_columns(
            pl.col("y").shift(1).over("mun").alias("y_l"), pl.col("wy").shift(1).over("mun").alias("wy_l"))
        P = (P.join(cov.select("mun", "dp"), on="mun", maintain_order="left").sort("mun", "anio").drop_nulls().filter(pl.col("y").is_not_nan() & pl.col("y_l").is_not_nan())
             .with_columns((pl.col("dp").cast(pl.Utf8) + "_" + pl.col("anio").cast(pl.Utf8)).alias("dpa")))
        out[c] = {nom: feols_dos_vias(P, "y", ["y_l", "wy_l"], fe, "dp") for nom, fe in (("anio", ("mun", "anio")), ("dpto_anio", ("mun", "dpa")))}
        print(f"  panel {c:>6} δ={out[c]['anio']['delta']} | dpto×año δ={out[c]['dpto_anio']['delta']}", flush=True)
    _json(out, rutas.OUT / "panel.json")


def etapa_departamental(nsim: int = 999):
    """I de Moran de la SMR departamental de 2022 (sin suavizar, como el trabajo previo) para las 14 series, con la
    armonización de edades de este repositorio y la contigüidad departamental derivada del MGN."""
    D = pl.read_parquet(rutas.INTER / "departamental_2022.parquet")
    dps = sorted(D["dp"].unique().to_list())
    Wd = pesos.matriz_departamentos(dps)
    rng = np.random.default_rng(SEMILLA)
    sel = pl.read_csv(rutas.CONFIG / "seleccion_causas_departamental.csv")
    previo = {str(r["codigo_667"]): r["I_moran_2022"] for r in sel.iter_rows(named=True)}
    out = {}
    for k in CLAVES:
        c = clave(k)
        d = pl.DataFrame({"dp": dps}).join(D.filter(pl.col("causa") == c), on="dp", how="left", maintain_order="left")
        smr = (d["O"] / d["E"]).to_numpy()
        I = moran(smr, Wd)
        sims = np.array([moran(rng.permutation(smr), Wd) for _ in range(nsim)])
        out[c] = {"I": I, "p_perm": float((np.sum(sims >= I) + 1) / (nsim + 1)), "I_trabajo_previo": previo.get(c)}
    _json({"meta": {"departamentos": dps, "nsim": nsim}, "causas": out}, rutas.OUT / "departamental.json")
    for c, v in out.items():
        print(f"  dpto {c:>6} I={v['I']:.3f} (previo {v['I_trabajo_previo']})")
