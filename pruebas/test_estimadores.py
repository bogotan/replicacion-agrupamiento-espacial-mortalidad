"""Estimadores escritos a mano (src/agrupamiento/estimadores.py) contra casos conocidos o contra paquetes."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
from agrupamiento import estimadores as est  # noqa: E402


def _reticula(k=15):
    """Contigüidad tipo torre en una retícula k × k: binaria simétrica y estandarizada por filas."""
    n = k * k; A = np.zeros((n, n))
    for i in range(k):
        for j in range(k):
            a = i * k + j
            for di, dj in ((1, 0), (0, 1)):
                if i + di < k and j + dj < k:
                    b = (i + di) * k + j + dj; A[a, b] = A[b, a] = 1
    return A, A / A.sum(1, keepdims=True)


def test_moran_contra_esda():
    import esda
    import libpysal
    A, Wd = _reticula(12)
    w = libpysal.weights.lat2W(12, 12, rook=True); w.transform = "r"
    x = np.random.default_rng(1).standard_normal(144).cumsum()
    assert est.moran(x, Wd) == pytest.approx(esda.Moran(x, w, permutations=0).I, abs=1e-12)
    assert est.moran_columnas(x[:, None], Wd)[0] == pytest.approx(est.moran(x, Wd), abs=1e-12)


def test_esperanza_moran_h0_solo_constante():
    _, Wd = _reticula(10)
    assert est.esperanza_moran_h0(np.empty((100, 0)), Wd) == pytest.approx(-1 / 99, abs=1e-12)


def test_shapley_suma_y_simetria():
    I = {"0": 0.6, "S": 0.5, "G": 0.45, "Q": 0.55, "SG": 0.35, "SQ": 0.42, "GQ": 0.40, "SGQ": 0.30}
    phi = est.shapley(I)
    assert sum(phi.values()) == pytest.approx(I["0"] - I["SGQ"], abs=1e-12)
    aditivo = {"0": 1.0, "S": 0.9, "G": 0.8, "Q": 0.7, "SG": 0.7, "SQ": 0.6, "GQ": 0.5, "SGQ": 0.4}
    assert est.shapley(aditivo) == pytest.approx({"S": 0.1, "G": 0.2, "Q": 0.3}, abs=1e-12)


def test_bayes_empirico():
    E = np.array([1.0, 5.0, 50.0, 500.0]); O = np.array([3.0, 4.0, 60.0, 480.0])
    th, info = est.bayes_empirico(O, E)
    r = O / E; m = O.sum() / E.sum(); s2 = (E * (r - m) ** 2).sum() / E.sum() - m / E.mean()
    C = s2 / (s2 + m / E)
    assert th == pytest.approx(m + C * (r - m)) and info["s2"] == pytest.approx(s2)
    assert abs(th[0] - m) < abs(r[0] - m) and abs(th[3] - r[3]) < abs(th[0] - r[0])     # más contracción si E es chico


def test_efectos_lesage_formula_espectral_contra_inversa():
    _, Wd = _reticula(10)
    ev = np.linalg.eigvals(Wd)
    beta, theta, rho = 0.4, 0.25, 0.55
    d_inv, i_inv, t_inv = est.efectos_sdm_puntuales(beta, theta, rho, Wd)
    d_esp = beta * np.real(np.mean(1 / (1 - rho * ev))) + theta * np.real(np.mean(ev / (1 - rho * ev)))
    assert d_esp == pytest.approx(d_inv, abs=1e-10)
    assert t_inv == pytest.approx((beta + theta) / (1 - rho), abs=1e-10)


def test_slx_errores_agrupados_contra_statsmodels():
    import statsmodels.api as sm
    rng = np.random.default_rng(3)
    _, Wd = _reticula(12); n = 144
    X = rng.standard_normal((n, 2)); g = np.repeat(np.arange(12), 12)
    y = 0.5 + X @ [0.3, -0.2] + Wd @ X @ [0.1, 0.2] + rng.standard_normal(n) + np.repeat(rng.standard_normal(12), 12)
    w = rng.uniform(0.5, 2, n)
    la = np.repeat(np.linspace(0, 1, 12), 12); lo = np.tile(np.linspace(0, 1, 12), 12)
    r = est.slx_inferencia(y, X, Wd, w, g, la, lo)
    Z = np.column_stack([np.ones(n), X, Wd @ X])
    m = sm.WLS(y, Z, weights=w).fit(cov_type="cluster", cov_kwds={"groups": g})
    assert r["dir"] == pytest.approx(m.params[1:3].tolist(), abs=1e-10)
    assert r["cl"]["se_dir"] == pytest.approx(m.bse[1:3].tolist(), rel=1e-8)
    V = m.cov_params(); se_t = [np.sqrt(V[1 + j, 1 + j] + V[3 + j, 3 + j] + 2 * V[1 + j, 3 + j]) for j in range(2)]
    assert r["cl"]["se_tot"] == pytest.approx(se_t, rel=1e-8)


def test_bh_contra_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = np.random.default_rng(5).uniform(0, 0.2, 50)
    assert est.bh(p) == pytest.approx(multipletests(p, method="fdr_bh")[1], abs=1e-12)


def test_feols_dos_vias_contra_dummies():
    import polars as pl
    rng = np.random.default_rng(7)
    u = np.repeat(np.arange(40), 6); t = np.tile(np.arange(6), 40); g = u // 5
    x1, x2 = rng.standard_normal(240), rng.standard_normal(240)
    y = 0.4 * x1 - 0.3 * x2 + rng.standard_normal(40)[u] + rng.standard_normal(6)[t] + rng.standard_normal(240)
    P = pl.DataFrame({"y": y, "x1": x1, "x2": x2, "u": u, "t": t, "g": g})
    r = est.feols_dos_vias(P, "y", ["x1", "x2"], ("u", "t"), "g")
    D = np.column_stack([x1, x2] + [(u == k).astype(float) for k in range(40)] + [(t == k).astype(float) for k in range(1, 6)])
    b = np.linalg.lstsq(D, y, rcond=None)[0]
    assert r["gamma"] == pytest.approx(b[0], abs=1e-8) and r["delta"] == pytest.approx(b[1], abs=1e-8)
    # y sin variación dentro de los efectos fijos (como la leucemia en el panel anual): no estimable, sin ruido
    P0 = P.with_columns(pl.Series("y", rng.standard_normal(40)[u] + rng.standard_normal(6)[t]))
    r0 = est.feols_dos_vias(P0, "y", ["x1", "x2"], ("u", "t"), "g")
    assert r0["estimable"] is False and r0["delta"] is None and r["estimable"] is True


def test_escala_icar_media_geometrica_uno():
    A, _ = _reticula(8)
    ev, U, c = est.escala_icar(A)
    Rmas = U[:, ev > 0] @ np.diag(1 / (c * ev[ev > 0])) @ U[:, ev > 0].T
    assert np.exp(np.mean(np.log(np.diag(Rmas)))) == pytest.approx(1.0, abs=1e-10)


def test_bym2_sin_efecto_espacial_coincide_con_poisson():
    """Sin variación extra-Poisson, los coeficientes del BYM2 coinciden con los de un GLM de Poisson (statsmodels)."""
    import statsmodels.api as sm
    rng = np.random.default_rng(11)
    A, _ = _reticula(12); n = 144
    ev, U, c = est.escala_icar(A)
    X = rng.standard_normal((n, 1)); E = rng.uniform(20, 60, n)
    O = rng.poisson(E * np.exp(-0.2 + 0.3 * X[:, 0])).astype(float)
    r = est.bym2_laplace(O, E, X, ev, U, c)
    glm = sm.GLM(O, sm.add_constant(X), family=sm.families.Poisson(), offset=np.log(E)).fit()
    assert r["sigma2"] < 0.01
    assert np.array(r["beta"]) == pytest.approx(glm.params, abs=0.02)


def test_bym2_recupera_parametros_simulados():
    """Caso conocido: efecto BYM2 con σ = 0,5 y φ = 0,8 en una retícula de 20 × 20; β se recupera dentro de tres
    errores estándar y φ queda por encima de 0,5."""
    rng = np.random.default_rng(13)
    A, _ = _reticula(20); n = 400
    ev, U, c = est.escala_icar(A)
    sig, phi = 0.5, 0.8
    lam_mas = np.where(ev > 1e-9, 1 / (c * np.where(ev > 1e-9, ev, 1)), 0.0)
    u_est = U @ (np.sqrt(lam_mas) * rng.standard_normal(n))
    b = sig * (np.sqrt(1 - phi) * rng.standard_normal(n) + np.sqrt(phi) * u_est)
    X = rng.standard_normal((n, 1)); E = rng.uniform(5, 30, n)
    O = rng.poisson(E * np.exp(0.1 + 0.25 * X[:, 0] + b)).astype(float)
    r = est.bym2_laplace(O, E, X, ev, U, c)
    se = np.sqrt(np.diag(np.array(r["cov_beta"])))
    assert abs(r["beta"][1] - 0.25) < 3 * se[1]
    assert r["phi"] > 0.5 and 0.1 < r["sigma2"] < 0.6
