"""Figuras (en color y en escala de grises) y tablas CSV del artículo, desde output/*.json y data/intermedios/.

Las figuras se diseñan para leerse en grises (la revista las exige así): las series se distinguen por la forma del
marcador y el relleno. Se guardan en output/figuras/ y output/figuras/gris/.
"""
from __future__ import annotations

import json
import math

import numpy as np
import polars as pl

from . import rutas
from .constantes import BLOQUES, COVS, EXCLUIR, NOMBRE_COV, clave
from .estimadores import bh as _bh
from .estimadores import bayes_empirico


def smr_eb() -> pl.DataFrame:
    cov = pl.read_parquet(rutas.INTER / "covariables.parquet").sort("mun")
    oe = pl.read_parquet(rutas.INTER / "observados_esperados.parquet")
    filas = []
    for c in oe["causa"].unique(maintain_order=True).to_list():
        d = cov.select("mun").join(oe.filter(pl.col("causa") == c), on="mun", how="left", maintain_order="left")
        th, _ = bayes_empirico(d["O"].to_numpy().astype(float), d["E"].to_numpy())
        filas.append(d.with_columns(pl.Series("smr_eb", th)))
    return pl.concat(filas)


AZUL, NAR, GRIS, VERDE, MOR = "#11607A", "#C0601E", "#8a8a8a", "#1BA578", "#7B4FA3"


def coma(x, k=2):
    from decimal import Decimal, ROUND_HALF_UP
    s = str(Decimal(str(float(x))).quantize(Decimal(1).scaleb(-k), rounding=ROUND_HALF_UP))
    ent, _, dec = s.partition(".")
    neg = ent.startswith("-"); ent = ent.lstrip("-")
    ent = f"{int(ent):,}".replace(",", ".")
    return ("−" if neg else "") + ent + ("," + dec if dec else "")


def estilo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fam = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({"font.family": fam, "font.size": 9.5, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#9a9a9a", "axes.grid": True, "grid.color": "#e3e3e3", "grid.linewidth": 0.7,
                         "axes.axisbelow": True, "savefig.dpi": 300, "savefig.bbox": "tight", "figure.facecolor": "white",
                         "axes.formatter.use_locale": False})
    return plt


def figuras_y_tablas():
    """Figuras pensadas para leerse en escala de grises (la revista las exige así): las series se distinguen por
    la forma del marcador y el relleno, no solo por el color. Se guardan en color y en grises (figuras/gris/)."""
    plt = estilo()
    import geopandas as gpd
    import textwrap
    from matplotlib.ticker import FuncFormatter
    from matplotlib.lines import Line2D
    from matplotlib.colors import LogNorm, LinearSegmentedColormap
    R = json.loads((rutas.OUT / "resultados.json").read_text())
    C = R["causas"]
    RV = json.loads((rutas.OUT / "revision.json").read_text())
    BY = json.loads((rutas.OUT / "bym2.json").read_text())["causas"]
    FIG = rutas.FIG; TAB = rutas.TAB
    FIG.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
    for f in FIG.glob("*.png"):
        f.unlink()
    fmt2 = FuncFormatter(lambda v, _: coma(v, 2))
    fmt0 = FuncFormatter(lambda v, _: coma(v, 0))
    orden = sorted([c for c in C if c != "todas"], key=lambda c: -C[c]["moran_y"]["I"]) + ["todas"]
    et = {c: C[c]["etiqueta"] for c in C}
    ys = np.arange(len(orden))[::-1]

    # Figura 1: mapas de la SMR-EB (escala secuencial logarítmica, legible en grises)
    g = gpd.read_file(rutas.EXT63 / "mgn2018_mpio_cnpv.geojson")
    g["mun"] = g["MPIO_CDPMP"].astype(int)
    se = smr_eb()
    cmap = LinearSegmentedColormap.from_list("sec", ["#f7f7f7", "#bcd3dc", AZUL, "#082a35"])
    norm = LogNorm(vmin=0.4, vmax=2.5)
    sel = ["201", "604", "302", "512"]
    fig, axs = plt.subplots(1, 4, figsize=(13, 4.6))
    for ax, c in zip(axs, sel):
        d = se.filter(pl.col("causa") == c).select("mun", "smr_eb").to_pandas()
        gg = g.merge(d, on="mun", how="left")
        gg = gg[~gg["mun"].isin(EXCLUIR)]
        gg.assign(v=gg["smr_eb"].clip(0.4, 2.5)).plot(column="v", ax=ax, cmap=cmap, norm=norm, linewidth=0.05, edgecolor="#9a9a9a")
        ax.set_axis_off()
        ax.set_title(f"{textwrap.fill(et[c], 34)}\nI de Moran = {coma(C[c]['moran_y']['I'], 2)}", fontsize=9.5)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axs, orientation="horizontal", fraction=0.04, pad=0.02, aspect=50,
                      ticks=[0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5], format=FuncFormatter(lambda v, _: coma(v, 1)))
    cb.set_label("Razón estandarizada de mortalidad suavizada (Bayes empírico), 2018–2025; 1 = promedio nacional; escala logarítmica")
    fig.savefig(FIG / "fig1_mapas_smr.png", facecolor="white"); plt.close(fig)

    # Figura 2: I inicial, I residual real, distribución placebo y efectos fijos
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    ei_fe = RV["EI_H0"]["SGQ+FE"]
    for yi, c in zip(ys, orden):
        d = RV["causas"][c]
        lo, med, hi = d["I_SGQ_placebo_p05_p50_p95"]
        ax.plot([lo, hi], [yi, yi], color="#9a9a9a", lw=3.2, alpha=0.55, solid_capstyle="butt", zorder=1)
        ax.scatter(d["I0"], yi, s=36, facecolor="white", edgecolor="#222222", lw=1.1, zorder=3, marker="o")
        ax.scatter(med, yi, s=34, color="#6b6b6b", marker="|", lw=1.6, zorder=3)
        ax.scatter(d["I_SGQ"], yi, s=36, color=AZUL, marker="o", zorder=4)
        ax.scatter(d["I_SGQFE"], yi, s=34, color=NAR, marker="D", zorder=4)
    ax.axvline(0, color="#6b6b6b", lw=0.8)
    ax.axvline(ei_fe, color="#6b6b6b", lw=0.8, ls=(0, (2, 2)))
    ax.set_yticks(ys); ax.set_yticklabels([et[c] for c in orden]); ax.xaxis.set_major_formatter(fmt2)
    ax.set_xlabel("I de Moran (contigüidad tipo reina)")
    hs = [Line2D([], [], marker="o", ls="", mfc="white", mec="#222222", label="log SMR-EB, sin covariables"),
          Line2D([], [], color="#9a9a9a", lw=3.2, alpha=0.55, label="residuo con 16 covariables placebo (P5–P95)"),
          Line2D([], [], marker="o", ls="", color=AZUL, label="residuo con S, G y Q"),
          Line2D([], [], marker="D", ls="", color=NAR, label="residuo con S, G, Q y efectos fijos"),
          Line2D([], [], color="#6b6b6b", lw=0.8, ls=(0, (2, 2)), label="E[I] bajo H0 con efectos fijos")]
    ax.legend(handles=hs, loc="lower right", frameon=False, fontsize=8)
    ax.grid(axis="y", visible=False)
    fig.savefig(FIG / "fig2_moran_placebo.png", facecolor="white"); plt.close(fig)

    # Figura 3: contribución de Shapley de cada bloque por encima del placebo de igual tamaño
    fig, axs = plt.subplots(1, 3, figsize=(12, 6.2), sharey=True)
    nomb = {"S": "Sociodemográfico (7 variables)", "G": "Geografía y oferta (5)", "Q": "Calidad del registro (4)"}
    for ax, b in zip(axs, "SGQ"):
        for yi, c in zip(ys, orden):
            d = RV["causas"][c]["bloques"][b]
            ex = 100 * d["shapley_exceso"]; lo, hi = (100 * x for x in d["shapley_exceso_ic"])
            fuera = lo > 0
            ax.plot([lo, hi], [yi, yi], color="#222222" if fuera else "#9a9a9a", lw=1.3)
            ax.scatter(ex, yi, s=26, marker="o", facecolor="#222222" if fuera else "white", edgecolor="#222222" if fuera else "#7a7a7a", zorder=3)
        ax.axvline(0, color="#6b6b6b", lw=0.8)
        ax.set_title(nomb[b], fontsize=9.5); ax.xaxis.set_major_formatter(fmt0); ax.grid(axis="y", visible=False)
        ax.set_xlim(-60, 40)
    axs[0].set_yticks(ys); axs[0].set_yticklabels([et[c] for c in orden])
    fig.supxlabel("Contribución de Shapley menos la de un bloque placebo del mismo tamaño (% de la I inicial), intervalo P2,5–P97,5 del placebo", fontsize=9)
    fig.savefig(FIG / "fig3_shapley_placebo.png", facecolor="white"); plt.close(fig)

    # Figura 4: BYM2, varianza del efecto aleatorio y fracción estructurada
    fig, axs = plt.subplots(1, 2, figsize=(11, 6.0), sharey=True, gridspec_kw={"width_ratios": [1.5, 1]})
    for yi, c in zip(ys, orden):
        b0, b1 = BY[c]["nulo"], BY[c]["completo"]
        axs[0].plot([b1["sigma2"], b0["sigma2"]], [yi, yi], color="#cfcfcf", lw=1.4)
        axs[0].scatter(b0["sigma2"], yi, s=30, facecolor="white", edgecolor="#222222", zorder=3)
        axs[0].scatter(b1["sigma2"], yi, s=30, color=AZUL, zorder=3)
        axs[1].plot([b1["phi"], b0["phi"]], [yi, yi], color="#cfcfcf", lw=1.4)
        axs[1].scatter(b0["phi"], yi, s=30, facecolor="white", edgecolor="#222222", zorder=3)
        axs[1].scatter(b1["phi"], yi, s=30, color=AZUL, zorder=3)
    axs[0].set_xscale("log"); axs[0].xaxis.set_major_formatter(FuncFormatter(lambda v, _: coma(v, 2) if v < 1 else coma(v, 0)))
    axs[0].set_xlabel("Varianza marginal σ² del efecto BYM2 (escala logarítmica)")
    axs[1].set_xlabel("Fracción estructurada φ"); axs[1].set_xlim(0, 1.02); axs[1].xaxis.set_major_formatter(FuncFormatter(lambda v, _: coma(v, 1)))
    axs[0].set_yticks(ys); axs[0].set_yticklabels([et[c] for c in orden])
    for ax in axs:
        ax.grid(axis="y", visible=False)
    axs[1].legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec="#222222", label="sin covariables"),
                           Line2D([], [], marker="o", ls="", color=AZUL, label="con S, G, Q y sus rezagos")],
                  frameon=False, loc="lower left", fontsize=8)
    fig.savefig(FIG / "fig4_bym2.png", facecolor="white"); plt.close(fig)

    # Figura 5: efecto total en las especificaciones coprincipales y en el SDEM
    filas = efectos_seleccionados(R, RV, BY)
    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    yy = np.arange(len(filas))[::-1]
    mk = {"sdem": ("o", "#9a9a9a", "SDEM (máxima verosimilitud, sin pesos)"), "mcp": ("s", AZUL, "SLX ponderado, errores agrupados por departamento"),
          "bym": ("^", NAR, "Poisson con BYM2 y rezagos de las covariables")}
    for yi, f in zip(yy, filas):
        for j, (kk, (m, col, _)) in enumerate(mk.items()):
            v, lo, hi = f[kk]
            off = (j - 1) * 0.22
            sig = lo > 0 or hi < 0
            ax.plot([lo, hi], [yi + off, yi + off], color=col, lw=1.2)
            ax.scatter(v, yi + off, marker=m, s=26, facecolor=col if sig else "white", edgecolor=col, zorder=3)
    ax.axvline(0, color="#6b6b6b", lw=0.8)
    ax.set_xscale("symlog", linthresh=50)
    ax.set_xticks([-50, -25, 0, 25, 50, 100, 200]); ax.xaxis.set_major_formatter(fmt0)
    ax.set_yticks(yy); ax.set_yticklabels([f["rotulo"] for f in filas], fontsize=8.6)
    ax.set_xlabel("Efecto total de una desviación estándar (PDET: pertenecer) sobre la SMR, %, IC del 95 %")
    ax.legend(handles=[Line2D([], [], marker=m, ls="", color=col, label=l) for m, col, l in mk.values()],
              frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.09), ncol=1, fontsize=8)
    ax.grid(axis="y", visible=False)
    fig.savefig(FIG / "fig5_efectos_coprincipales.png", facecolor="white"); plt.close(fig)

    from PIL import Image
    (FIG / "gris").mkdir(exist_ok=True)
    for f in (FIG / "gris").glob("*.png"):
        f.unlink()
    for f in FIG.glob("*.png"):
        Image.open(f).convert("L").save(FIG / "gris" / f.name)
    escribir_tablas(R, TAB)
    escribir_tablas_revision(R, RV, BY, TAB)
    print("  -> figuras y tablas")


PARES = [("201", "altitud", "Cáncer de estómago · altitud"), ("103", "altitud", "Vectores · altitud"),
         ("404", "altitud", "Dificultad resp. del RN · altitud"), ("404", "nbi", "Dificultad resp. del RN · NBI"),
         ("108", "indigena", "Neumonía e IRA · indígena"), ("204", "internet", "Cáncer de pulmón · internet"),
         ("204", "superior", "Cáncer de pulmón · educación superior"), ("210", "afro", "Cáncer de próstata · afrocolombiana"),
         ("604", "internet", "Alzheimer y otras · internet"), ("604", "completitud", "Alzheimer y otras · completitud"),
         ("todas", "sin_asistencia", "Todas las causas · sin asistencia"), ("todas", "basura", "Todas las causas · códigos basura"),
         ("302", "basura", "Hipertensivas · códigos basura"), ("302", "basura_sin_hta", "Hipertensivas · basura sin I10–I15"),
         ("512", "superior", "Homicidios · educación superior"), ("512", "pdet", "Homicidios · PDET")]


def efectos_seleccionados(R, RV, BY):
    """Efecto total (β+θ) en % con IC 95 %: SDEM-ML, SLX-MCP con errores agrupados y BYM2 con WX."""
    k = len(COVS)
    filas = []
    for c, v, rot in PARES:
        sin_hta = v == "basura_sin_hta"
        j = COVS.index("basura" if sin_hta else v)
        f = {"causa": c, "variable": v, "rotulo": rot}
        if sin_hta:
            f["sdem"] = (np.nan, np.nan, np.nan)
        else:
            im = R["causas"][c]["impactos"]["SDEM"]["total"]
            f["sdem"] = tuple(100 * np.expm1(np.array([im["media"][j], im["lo"][j], im["hi"][j]])))
        m = RV["slx"][c]["MCP_sin_hta" if sin_hta else "MCP"]
        t, s = m["tot"][j], m["cl"]["se_tot"][j]
        f["mcp"] = tuple(100 * np.expm1(np.array([t, t - 1.96 * s, t + 1.96 * s])))
        f["mcp_q"] = m["cl"].get("q_tot", [np.nan] * k)[j]
        bm = BY[c]["sin_hta"] if sin_hta else BY[c]["completo"]
        bb = np.array(bm["beta"]); V = np.array(bm["cov_beta"])
        tb = bb[1 + j] + bb[1 + k + j]; sb = math.sqrt(V[1 + j, 1 + j] + V[1 + k + j, 1 + k + j] + 2 * V[1 + j, 1 + k + j])
        f["bym"] = tuple(100 * np.expm1(np.array([tb, tb - 1.96 * sb, tb + 1.96 * sb])))
        db = bb[1 + j]; sd = math.sqrt(V[1 + j, 1 + j])
        f["bym_dir"] = tuple(100 * np.expm1(np.array([db, db - 1.96 * sd, db + 1.96 * sd])))
        filas.append(f)
    return filas


def escribir_tablas_revision(R, RV, BY, TAB):
    from scipy.stats import norm
    k = len(COVS)
    # tabla 7: placebo, numerador y denominador
    filas = []
    for c, d in RV["causas"].items():
        filas.append({"causa": c, "etiqueta": R["causas"][c]["etiqueta"], "I0": d["I0"], "I_SGQ": d["I_SGQ"], "I_SGQ_FE": d["I_SGQFE"],
                      "caida_real": d["caida_real"], "caida_placebo_media": d["caida_placebo_media"], "caida_placebo_p95": d["caida_placebo_p95"],
                      "caida_exceso": d["caida_exceso"], "p_placebo": d["p_placebo"],
                      "I_SGQ_placebo_p05": d["I_SGQ_placebo_p05_p50_p95"][0], "I_SGQ_placebo_p95": d["I_SGQ_placebo_p05_p50_p95"][2],
                      "caida_numerador": d["caida_num"], "caida_denominador": d["caida_den"],
                      **{f"sh_{b}_real": d["bloques"][b]["shapley_real"] for b in "SGQ"},
                      **{f"sh_{b}_placebo": d["bloques"][b]["shapley_placebo_media"] for b in "SGQ"},
                      **{f"sh_{b}_exceso": d["bloques"][b]["shapley_exceso"] for b in "SGQ"},
                      **{f"sh_{b}_lo": d["bloques"][b]["shapley_exceso_ic"][0] for b in "SGQ"},
                      **{f"sh_{b}_hi": d["bloques"][b]["shapley_exceso_ic"][1] for b in "SGQ"},
                      **{f"sh_{b}_p": d["bloques"][b]["p_placebo"] for b in "SGQ"},
                      "pares_no_monotonos": len(d["no_monotono"])})
    pl.DataFrame(filas).write_csv(TAB / "tabla7_placebo_shapley.csv")
    # tabla 8: efectos totales en las tres especificaciones, todas las combinaciones
    filas = []
    pb = []
    for c in R["causas"]:
        bm = BY[c]["completo"]; bb = np.array(bm["beta"]); V = np.array(bm["cov_beta"])
        for j, v in enumerate(COVS):
            im = R["causas"][c]["impactos"]["SDEM"]["total"]
            m = RV["slx"][c]["MCP"]; mo = RV["slx"][c]["MCO"]
            tb = bb[1 + j] + bb[1 + k + j]; sb = math.sqrt(V[1 + j, 1 + j] + V[1 + k + j, 1 + k + j] + 2 * V[1 + j, 1 + k + j])
            pbym = 2 * (1 - norm.cdf(abs(tb / sb))); pb.append(pbym)
            filas.append({"causa": c, "variable": v, "sdem_total": im["media"][j], "sdem_lo": im["lo"][j], "sdem_hi": im["hi"][j],
                          "mco_total": mo["tot"][j], "mco_se_cl": mo["cl"]["se_tot"][j], "mco_q_cl": mo["cl"]["q_tot"][j],
                          "mcp_directo": m["dir"][j], "mcp_se_dir_cl": m["cl"]["se_dir"][j],
                          "mcp_total": m["tot"][j], "mcp_se_cl": m["cl"]["se_tot"][j], "mcp_se_conley": m["conley"]["se_tot"][j],
                          "mcp_q_cl": m["cl"]["q_tot"][j], "mcp_q_conley": m["conley"]["q_tot"][j],
                          "bym_total": tb, "bym_se": sb, "bym_directo": bb[1 + j], "bym_se_dir": math.sqrt(V[1 + j, 1 + j]), "bym_p": pbym})
    T = pl.DataFrame(filas).with_columns(pl.Series("bym_q", _bh(pb)))
    T = T.with_columns(((pl.col("mcp_q_cl") < 0.05) & (pl.col("bym_q") < 0.05) & (pl.col("mcp_total").sign() == pl.col("bym_total").sign())
                        & (pl.col("sdem_total").sign() == pl.col("mcp_total").sign())).alias("sobrevive"))
    T.write_csv(TAB / "tabla8_efectos_coprincipales.csv")
    # tabla 9: BYM2
    filas = []
    for c, b in BY.items():
        filas.append({"causa": c, "etiqueta": R["causas"][c]["etiqueta"], "sigma2_nulo": b["nulo"]["sigma2"], "phi_nulo": b["nulo"]["phi"],
                      "frontera_nulo": b["nulo"]["frontera"], "sigma2_completo": b["completo"]["sigma2"], "phi_completo": b["completo"]["phi"],
                      "frontera_completo": b["completo"]["frontera"], "caida_sigma2": 1 - b["completo"]["sigma2"] / b["nulo"]["sigma2"],
                      "moran_u_nulo": b["moran_u_nulo"], "moran_u_completo": b["moran_u_completo"],
                      "conv_nulo": b["nulo"]["conv"], "conv_completo": b["completo"]["conv"]})
    pl.DataFrame(filas).write_csv(TAB / "tabla9_bym2.csv")
    # tabla 10: panel anual
    if (rutas.OUT / "panel.json").exists():
        P = json.loads((rutas.OUT / "panel.json").read_text())
        pl.DataFrame([{"causa": c, "etiqueta": R["causas"][c]["etiqueta"], **{f"{m}_{kk}": v for m, d in r.items() for kk, v in d.items()}}
                      for c, r in P.items()]).write_csv(TAB / "tabla10_panel.csv")




def elegido(r):
    """Modelo con menor AIC entre SDM y SDEM (los dos anidan SAR/SEM/SLX)."""
    return "SDM" if r["modelos"]["SDM"]["aic"] < r["modelos"]["SDEM"]["aic"] else "SDEM"


def escribir_tablas(R, TAB):
    C = R["causas"]
    cov = pl.read_parquet(rutas.INTER / "covariables.parquet")
    # Tabla 1: descriptivos de covariables
    filas = []
    for v in COVS:
        x = cov[v].drop_nulls().to_numpy().astype(float)
        if v in ("ldens", "ldist_capital", "ldist_uci", "lcamas"):
            base = {"ldens": None, "ldist_capital": "dist_capital", "ldist_uci": "dist_uci", "lcamas": None}[v]
        filas.append({"variable": v, "nombre": NOMBRE_COV[v], "bloque": next(b for b, vs in BLOQUES.items() if v in vs),
                      "media": float(x.mean()), "de": float(x.std()), "p10": float(np.percentile(x, 10)),
                      "mediana": float(np.median(x)), "p90": float(np.percentile(x, 90)), "faltantes": int(cov[v].is_null().sum())})
    pl.DataFrame(filas).write_csv(TAB / "tabla1_covariables.csv")
    # Tabla de resultados por causa
    filas = []
    for c, r in C.items():
        m = elegido(r)
        filas.append({"causa": c, "etiqueta": r["etiqueta"], "defunciones": r["O"], "ceros": r["ceros"],
                      "I_y": r["moran_y"]["I"], "p_y": r["moran_y"]["p_sim"], "I_AR": r["moran_ar"]["I"],
                      "I_S": r["I_resid"]["S"]["reina"], "I_SG": r["I_resid"]["SG"]["reina"], "I_SGQ": r["I_resid"]["SGQ"]["reina"],
                      "z_SGQ": r["ols"]["moran_res"][1], "p_SGQ": r["ols"]["moran_res"][2],
                      "I_FE": r["I_resid"]["FE"]["reina"], "I_SGQ_FE": r["I_resid"]["SGQ+FE"]["reina"],
                      "p_SGQ_FE": r["ols_fe"]["moran_res"][2],
                      "LMe": r["ols"]["lm_error"][0], "LMl": r["ols"]["lm_lag"][0], "RLMe": r["ols"]["rlm_error"][0],
                      "p_RLMe": r["ols"]["rlm_error"][1], "RLMl": r["ols"]["rlm_lag"][0], "p_RLMl": r["ols"]["rlm_lag"][1],
                      "aic_ols": r["modelos"]["OLS"]["aic"], "aic_sar": r["modelos"]["SAR"]["aic"], "aic_sem": r["modelos"]["SEM"]["aic"],
                      "aic_sdm": r["modelos"]["SDM"]["aic"], "aic_sdem": r["modelos"]["SDEM"]["aic"], "aic_slx": r["modelos"]["SLX"]["aic"],
                      "elegido": m, "rho_sdm": r["modelos"]["SDM"]["rho"], "lam_sdem": r["modelos"]["SDEM"]["lam"],
                      "lr_sdm_sar_p": r["lr"]["SDM_vs_SAR"][1], "lr_sdm_sem_p": r["lr"]["SDM_vs_SEM"][1], "lr_sdem_sem_p": r["lr"]["SDEM_vs_SEM"][1],
                      "lr_sdem_slx_p": r["lr"]["SDEM_vs_SLX"][1],
                      "sh_S": r["shapley"]["reina"]["S"], "sh_G": r["shapley"]["reina"]["G"], "sh_Q": r["shapley"]["reina"]["Q"],
                      "ar2_ols": r["ols"]["ar2"], "ar2_fe": r["ols_fe"]["ar2"],
                      "I_sinpand": r["sin_pandemia"]["moran_y"], "Ires_sinpand": r["sin_pandemia"]["moran_res"],
                      "chow_p": r.get("regimenes", {}).get("chow_global", [None, None])[1],
                      "rho_sdm_fe": r.get("fe_espacial", {}).get("rho_sdm"), "se_rho_sdm_fe": r.get("fe_espacial", {}).get("se_rho_sdm"),
                      "lam_sdem_fe": r.get("fe_espacial", {}).get("lam_sdem"), "se_lam_sdem_fe": r.get("fe_espacial", {}).get("se_lam_sdem"),
                      "aic_sdm_fe": r.get("fe_espacial", {}).get("aic_sdm"), "aic_sdem_fe": r.get("fe_espacial", {}).get("aic_sdem"),
                      "RLMe_fe": r["ols_fe"]["rlm_error"][0], "p_RLMe_fe": r["ols_fe"]["rlm_error"][1],
                      "RLMl_fe": r["ols_fe"]["rlm_lag"][0], "p_RLMl_fe": r["ols_fe"]["rlm_lag"][1],
                      "moran_res_sdem": r["moran_res_sdem"], "moran_res_sdm": r["moran_res_sdm"]})
    pl.DataFrame(filas).write_csv(TAB / "tabla2_moran_modelos.csv")
    # Impactos del SDEM (especificación principal) y del SDM, todas las causas y covariables, con q de Benjamini-Hochberg
    filas = []
    for m in ("SDEM", "SDM"):
        for c, r in C.items():
            for j, v in enumerate(COVS):
                d = {"modelo": m, "causa": c, "etiqueta": r["etiqueta"], "variable": v, "nombre": NOMBRE_COV[v],
                     "aic_menor": elegido(r) == m}
                for t in ("directo", "indirecto", "total"):
                    im = r["impactos"][m][t]
                    d[f"{t}"] = im["media"][j]; d[f"{t}_lo"] = im["lo"][j]; d[f"{t}_hi"] = im["hi"][j]; d[f"{t}_p"] = im["p"][j]
                filas.append(d)
    T = pl.DataFrame(filas)
    qs = []
    for m in ("SDEM", "SDM"):
        for t in ("directo", "total"):
            p = T.filter(pl.col("modelo") == m)[f"{t}_p"].to_numpy().astype(float).clip(1e-12, 1)
            m_ = len(p); o = np.argsort(p); q = np.empty(m_)
            q[o] = np.minimum.accumulate((p[o] * m_ / np.arange(1, m_ + 1))[::-1])[::-1].clip(max=1)
            qs.append((m, t, q))
    for t in ("directo", "total"):
        T = T.with_columns(pl.concat([pl.Series(q) for m, tt, q in qs if tt == t]).alias(f"{t}_q"))
    T.write_csv(TAB / "tabla3_impactos.csv")
    # Robustez por matriz de pesos
    filas = []
    for c, r in C.items():
        d = {"causa": c, "etiqueta": r["etiqueta"], "I_y_reina": r["moran_y"]["I"], "I_SGQ_reina": r["I_resid"]["SGQ"]["reina"]}
        for wn, o in r["otras_W"].items():
            d[f"I_y_{wn}"] = o["moran_y"]; d[f"I_SGQ_{wn}"] = o["moran_res"][0]; d[f"p_SGQ_{wn}"] = o["moran_res"][2]
            d[f"lam_sdem_{wn}"] = o["lam_sdem"]; d[f"rho_sdm_{wn}"] = o["rho_sdm"]
        filas.append(d)
    pl.DataFrame(filas).write_csv(TAB / "tabla4_robustez_W.csv")
    # MCO: coeficientes con EE clásicos, HC1 y Conley
    filas = []
    for c, r in C.items():
        for j, v in enumerate(COVS):
            filas.append({"causa": c, "variable": v, "beta": r["ols"]["beta"][j + 1], "se_ols": r["ols"]["se"][j + 1],
                          "se_hc1": r["conley"]["hc1"][j + 1], "se_conley50": r["conley"]["50"][j + 1],
                          "se_conley100": r["conley"]["100"][j + 1], "se_conley200": r["conley"]["200"][j + 1],
                          "beta_fe": (r["ols_fe"]["beta"][1 + r["ols_fe"]["nombres"].index(v)] if v in r["ols_fe"]["nombres"] else None),
                          "se_fe": (r["ols_fe"]["se"][1 + r["ols_fe"]["nombres"].index(v)] if v in r["ols_fe"]["nombres"] else None),
                          "beta_gmm_sdem": r["gmm_sdem"].get("beta", [None] * (j + 2))[j + 1],
                          "se_gmm_sdem": r["gmm_sdem"].get("se", [None] * (j + 2))[j + 1]})
    pl.DataFrame(filas).write_csv(TAB / "tabla6_mco_conley.csv")


