"""Estimadores escritos a mano (sin paquete) o envoltorios mínimos, con su ecuación del manuscrito.

  bayes_empirico()     ecuación (2): Bayes empírico global (Anselin, 2018)
  moran()              ecuación (3): I de Moran de un vector o de residuos
  esperanza_moran_h0() valor esperado de la I de los residuos bajo independencia, tr(MW)/(n−k)
  shapley()            ecuación (4): valor de Shapley de la caída de la I entre bloques
  impactos_sdm()       ecuación de los efectos de LeSage (directo, indirecto, total) por simulación
  impactos_sdem()      efectos del SDEM: β, θ y β+θ con errores analíticos
  conley_se(), conley_V(), cluster_se()  errores estándar de Conley (núcleo de Bartlett) y agrupados (CR1)
  slx_inferencia()     SLX por mínimos cuadrados (ponderados) con errores agrupados y de Conley
  bh()                 Benjamini-Hochberg
  feols_dos_vias()     ecuación (8): MCO con dos conjuntos de efectos fijos (proyecciones alternadas), errores agrupados
  bym2_laplace()       ecuación (9): Poisson con efecto aleatorio BYM2 escalado, por aproximación de Laplace
"""
from __future__ import annotations

import itertools
import math

import numpy as np

from .datos import hav


def bayes_empirico(O: np.ndarray, E: np.ndarray) -> tuple[np.ndarray, dict]:
    """Ecuación (2). Contracción de la SMR hacia la media global; s² por momentos, truncada en cero."""
    r = O / E
    m = O.sum() / E.sum()
    s2 = max(float((E * (r - m) ** 2).sum() / E.sum() - m / E.mean()), 0.0)
    C = s2 / (s2 + m / E)
    return m + C * (r - m), {"m": float(m), "s2": s2, "C_mediana": float(np.median(C)), "C": C}


def moran(e: np.ndarray, Wd: np.ndarray) -> float:
    """Ecuación (3) con W densa estandarizada por filas: I = n/S0 · e'We / e'e (e centrado)."""
    n = len(e); e = e - e.mean()
    return float((n / Wd.sum()) * (e @ Wd @ e) / (e @ e))


def moran_columnas(E_: np.ndarray, Wd: np.ndarray) -> np.ndarray:
    n = E_.shape[0]; E_ = E_ - E_.mean(0)
    return (n / Wd.sum()) * np.einsum("ij,ij->j", E_, Wd @ E_) / np.einsum("ij,ij->j", E_, E_)


def numerador_denominador(E_: np.ndarray, Wd: np.ndarray):
    n = E_.shape[0]; E_ = E_ - E_.mean(0)
    return np.einsum("ij,ij->j", E_, Wd @ E_) / n, np.einsum("ij,ij->j", E_, E_) / n


def esperanza_moran_h0(Xm: np.ndarray, Wd: np.ndarray) -> float:
    n = Wd.shape[0]
    Z = np.column_stack([np.ones(n), Xm]) if Xm.size else np.ones((n, 1))
    M = np.eye(n) - Z @ np.linalg.solve(Z.T @ Z, Z.T)
    return float(np.trace(M @ Wd) / (n - Z.shape[1]))


def zstd(v):
    v = np.asarray(v, float); return (v - v.mean()) / v.std()


def ols_resid(y, X):
    Z = np.column_stack([np.ones(len(y)), X]) if X.size else np.ones((len(y), 1))
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return y - Z @ b, b


def shapley(I: dict[str, float], bloques=("S", "G", "Q")) -> dict[str, float]:
    """Ecuación (4). I: I de los residuos por subconjunto de bloques, con claves «0», «S», «SG», …, «SGQ»
    (letras en el orden de `bloques`). Devuelve φ_k; Σ φ_k = I(∅) − I(todos)."""
    bl = list(bloques); m = len(bl); out = {}
    for b_ in bl:
        otros = [x for x in bl if x != b_]; s = 0.0
        for r_ in range(m):
            for sub in itertools.combinations(otros, r_):
                peso = math.factorial(len(sub)) * math.factorial(m - len(sub) - 1) / math.factorial(m)
                s_sin = "".join(x for x in bl if x in sub) or "0"
                s_con = "".join(x for x in bl if x in sub or x == b_)
                s += peso * (I[s_sin] - I[s_con])
        out[b_] = s
    return out


def conley_se(y, X, la, lo, cortes=(50, 100, 200)):
    Z = np.column_stack([np.ones(len(y)), X]); n, k = Z.shape
    b = np.linalg.solve(Z.T @ Z, Z.T @ y); e = y - Z @ b
    D = hav(la[:, None], lo[:, None], la[None, :], lo[None, :])
    iZZ = np.linalg.inv(Z.T @ Z)
    out = {}
    Ze = Z * e[:, None]
    for c in cortes:
        K = np.clip(1 - D / c, 0, None)
        out[c] = np.sqrt(np.diag(iZZ @ (Ze.T @ K @ Ze) @ iZZ * n / (n - k)))
    out["hc1"] = np.sqrt(np.diag(iZZ @ (Ze.T @ Ze) @ iZZ * n / (n - k)))
    return b, out


def impactos_sdm(m, k, eig, nsim, rng):
    """Efectos de LeSage para el SDM y el SAR: directo = (1/n) tr[(I−ρW)⁻¹(Iβ + Wθ)] = (1/n) Σ_j (β + θω_j)/(1 − ρω_j),
    total = (β+θ)/(1−ρ), indirecto = total − directo; inferencia con extracciones de la distribución asintótica."""
    b = np.asarray(m.betas).ravel(); V = np.asarray(m.vm)
    draws = rng.multivariate_normal(b, V, size=nsim, method="eigh")
    rho = draws[:, -1]
    beta = draws[:, 1:1 + k]
    theta = draws[:, 1 + k:1 + 2 * k] if draws.shape[1] > 2 + k else np.zeros_like(beta)
    s1 = np.array([np.real(np.mean(1 / (1 - r * eig))) for r in rho])
    s2 = np.array([np.real(np.mean(eig / (1 - r * eig))) for r in rho])
    dir_ = beta * s1[:, None] + theta * s2[:, None]
    tot = (beta + theta) / (1 - rho)[:, None]
    ind = tot - dir_
    res = {}
    for nom, arr in (("directo", dir_), ("indirecto", ind), ("total", tot)):
        res[nom] = {"media": arr.mean(0), "lo": np.percentile(arr, 2.5, 0), "hi": np.percentile(arr, 97.5, 0),
                    "p": 2 * np.minimum((arr > 0).mean(0), (arr < 0).mean(0))}
    return res


def efectos_sdm_puntuales(beta, theta, rho, Wd):
    """Directo, indirecto y total en el punto, por la inversa explícita (para probar impactos_sdm)."""
    n = Wd.shape[0]
    S = np.linalg.solve(np.eye(n) - rho * Wd, np.eye(n) * beta + Wd * theta)
    d = np.trace(S) / n; t = S.sum() / n
    return d, t - d, t


def impactos_sdem(m, k):
    from scipy.stats import norm
    b = np.asarray(m.betas).ravel(); V = np.asarray(m.vm)
    res = {}
    for nom, ii in {"directo": list(range(1, 1 + k)), "indirecto": list(range(1 + k, 1 + 2 * k))}.items():
        se = np.sqrt(np.diag(V)[ii]); mu = b[ii]
        res[nom] = {"media": mu, "lo": mu - 1.96 * se, "hi": mu + 1.96 * se, "p": 2 * (1 - norm.cdf(np.abs(mu / se)))}
    it = list(range(1, 1 + k)); jt = list(range(1 + k, 1 + 2 * k))
    mu = b[it] + b[jt]
    se = np.sqrt(np.array([V[i, i] + V[j, j] + 2 * V[i, j] for i, j in zip(it, jt)]))
    res["total"] = {"media": mu, "lo": mu - 1.96 * se, "hi": mu + 1.96 * se, "p": 2 * (1 - norm.cdf(np.abs(mu / se)))}
    return res


def bh(p):
    p = np.asarray(p, float).clip(1e-300, 1); m_ = len(p); o = np.argsort(p); q = np.empty(m_)
    q[o] = np.minimum.accumulate((p[o] * m_ / np.arange(1, m_ + 1))[::-1])[::-1].clip(max=1)
    return q


def cluster_se(Z, e, g):
    """Matriz de covarianza agrupada CR1 para MCO o MCP (Z y e ya multiplicados por la raíz del peso)."""
    n, k = Z.shape
    iZZ = np.linalg.inv(Z.T @ Z)
    M = np.zeros((k, k))
    for gg in np.unique(g):
        s = (Z[g == gg] * e[g == gg, None]).sum(0)
        M += np.outer(s, s)
    G = len(np.unique(g))
    return iZZ @ M @ iZZ * G / (G - 1) * (n - 1) / (n - k)


def conley_V(Z, e, la, lo, corte=100.0):
    n, k = Z.shape
    D = hav(la[:, None], lo[:, None], la[None, :], lo[None, :])
    K = np.clip(1 - D / corte, 0, None)
    iZZ = np.linalg.inv(Z.T @ Z); Ze = Z * e[:, None]
    return iZZ @ (Ze.T @ K @ Ze) @ iZZ * n / (n - k)


def slx_inferencia(y, X, Wd, w, g, la, lo):
    """SLX y = α + Xβ + WXθ + e por mínimos cuadrados con pesos w; errores agrupados por g y de Conley (100 km).
    Devuelve directo (β), indirecto (θ) y total (β+θ) con sus errores estándar."""
    Z = np.column_stack([np.ones(len(y)), X, Wd @ X]); k = X.shape[1]
    sw = np.sqrt(w); Zw = Z * sw[:, None]; yw = y * sw
    b = np.linalg.lstsq(Zw, yw, rcond=None)[0]; e = yw - Zw @ b
    out = {}
    for nom, V in (("cl", cluster_se(Zw, e, g)), ("conley", conley_V(Zw, e, la, lo))):
        se_d = np.sqrt(np.diag(V)[1:1 + k]); se_i = np.sqrt(np.diag(V)[1 + k:1 + 2 * k])
        se_t = np.sqrt(np.array([V[1 + j, 1 + j] + V[1 + k + j, 1 + k + j] + 2 * V[1 + j, 1 + k + j] for j in range(k)]))
        out[nom] = {"se_dir": se_d.tolist(), "se_ind": se_i.tolist(), "se_tot": se_t.tolist()}
    out["dir"] = b[1:1 + k].tolist(); out["ind"] = b[1 + k:1 + 2 * k].tolist(); out["tot"] = (b[1:1 + k] + b[1 + k:1 + 2 * k]).tolist()
    return out


def feols_dos_vias(P, yv, xv, fes, clus):
    """Ecuación (8): MCO con dos conjuntos de efectos fijos (proyecciones alternadas hasta 1e−10) y errores agrupados CR1."""
    M = P.select([yv] + xv).to_numpy().astype(float)
    grupos = [P[f].to_numpy() for f in fes]
    idx = [np.unique(gv, return_inverse=True)[1] for gv in grupos]
    D = M.copy()
    for _ in range(500):
        prev = D.copy()
        for ix in idx:
            sums = np.zeros((ix.max() + 1, D.shape[1])); np.add.at(sums, ix, D)
            D = D - (sums / np.bincount(ix)[:, None])[ix]
        if np.max(np.abs(D - prev)) < 1e-10:
            break
    y, X = D[:, 0], D[:, 1:]
    n = len(y)
    # Sin variación dentro de los efectos fijos (p. ej., el suavizado anual contrae todas las SMR al promedio), y queda
    # en ruido de redondeo y cualquier coeficiente sería ese ruido: no estimable.
    escala = np.linalg.norm(M - M.mean(0), axis=0)
    if np.any((escala < 1e-12) | (np.linalg.norm(D, axis=0) <= 1e-9 * escala)):
        return {"estimable": False, "gamma": None, "se_gamma": None, "delta": None, "se_delta": None, "n": int(n)}
    b = np.linalg.lstsq(X, y, rcond=None)[0]; e = y - X @ b
    k = X.shape[1] + sum(len(np.unique(gv)) for gv in grupos) - 1
    g = P[clus].to_numpy(); iXX = np.linalg.inv(X.T @ X); Mm = np.zeros((X.shape[1],) * 2)
    for gg in np.unique(g):
        sc = (X[g == gg] * e[g == gg, None]).sum(0); Mm += np.outer(sc, sc)
    G = len(np.unique(g)); V = iXX @ Mm @ iXX * G / (G - 1) * (n - 1) / (n - k)
    se = np.sqrt(np.diag(V))
    return {"estimable": True, "gamma": float(b[0]), "se_gamma": float(se[0]), "delta": float(b[1]), "se_delta": float(se[1]),
            "n": int(n)}


def escala_icar(A: np.ndarray):
    """Valores y vectores propios de R = D − A y factor de escala c tal que la media geométrica de diag((cR)⁺) sea 1."""
    R = np.diag(A.sum(1)) - A
    ev, U = np.linalg.eigh(R)
    ev = np.where(ev < 1e-9, 0.0, ev)
    diag_mas = (U[:, ev > 0] ** 2 / ev[ev > 0]).sum(1)
    return ev, U, float(np.exp(np.mean(np.log(diag_mas))))


def bym2_laplace(O, E, X, ev, U, c_esc, theta0=None):
    """Ecuación (9). Poisson con efecto aleatorio BYM2 (Riebler et al., 2016): b = σ(√(1−φ) v + √φ u*), u* ICAR escalado.
    Cov(b) = σ²[(1−φ)I + φ Q*⁺]; con R = UΛU', sus valores propios son σ²[(1−φ) + φ/(c·λ_j)] (λ_j > 0) y σ²(1−φ)
    en el modo constante. Hiperparámetros (log σ, logit φ) por máxima verosimilitud marginal (Laplace); β y b en la moda
    condicional; previa vaga N(0, 10²) sobre β. Varios arranques de Nelder-Mead y un reinicio desde el mejor."""
    from scipy.linalg import cho_factor, cho_solve
    from scipy.optimize import minimize
    from scipy.special import gammaln
    n = len(O); off = np.log(E)
    Z = np.column_stack([np.ones(n), X]) if X.size else np.ones((n, 1)); p = Z.shape[1]
    nulo = ev < 1e-9
    estado = {"z": np.concatenate([[np.log(O.sum() / E.sum())], np.zeros(p - 1 + n)])}
    lg = gammaln(O + 1).sum()

    def dvals(sig, phi):
        d = np.where(nulo, 1 - phi, (1 - phi) + phi / (c_esc * np.where(nulo, 1, ev)))
        return sig ** 2 * np.maximum(d, 1e-12)

    def obj(z, P):
        b, u = z[:p], z[p:]; eta = off + Z @ b + u
        return float(np.sum(O * eta - np.exp(eta)) - 0.5 * u @ P @ u - 0.5 * b @ b / 100.0)

    def hess(z, P):
        b, u = z[:p], z[p:]; mu = np.exp(off + Z @ b + u)
        H = np.empty((p + n, p + n))
        H[:p, :p] = Z.T @ (Z * mu[:, None]) + np.eye(p) / 100.0
        H[:p, p:] = (Z * mu[:, None]).T; H[p:, :p] = H[:p, p:].T; H[p:, p:] = P + np.diag(mu)
        return H, mu

    def modo(P):
        z = estado["z"].copy()
        for _ in range(50):
            H, mu = hess(z, P)
            b, u = z[:p], z[p:]
            gr = np.concatenate([Z.T @ (O - mu) - b / 100.0, (O - mu) - P @ u])
            paso = cho_solve(cho_factor(H), gr)
            t = 1.0; f0 = obj(z, P)
            while t > 1e-4 and obj(z + t * paso, P) < f0 - 1e-10:
                t /= 2
            z = z + t * paso
            if np.max(np.abs(t * paso)) < 1e-7:
                break
        H, _ = hess(z, P)
        return z, H

    def nlml(th):
        sig = np.exp(th[0]); phi = min(max(1 / (1 + np.exp(-th[1])), 1e-4), 0.999)
        d = dvals(sig, phi); P = (U / d) @ U.T
        z, H = modo(P); estado["z"] = z
        ld = 2 * np.log(np.diag(np.linalg.cholesky(H))).sum()
        return -(obj(z, P) - lg - 0.5 * np.log(d).sum() - 0.5 * p * np.log(100.0) - 0.5 * ld)

    opts = {"xatol": 1e-3, "fatol": 1e-3, "maxiter": 200}
    inicios = [np.array([np.log(0.4), 0.0]), np.array([np.log(0.4), 2.5])]
    if theta0 is not None:
        inicios.insert(0, np.asarray(theta0, float))
    r = min((minimize(nlml, t0, method="Nelder-Mead", options=opts) for t0 in inicios), key=lambda x: x.fun)
    r = minimize(nlml, r.x, method="Nelder-Mead", options=opts)
    sig = float(np.exp(r.x[0])); phi = float(min(max(1 / (1 + np.exp(-r.x[1])), 1e-4), 0.999))
    d = dvals(sig, phi); P = (U / d) @ U.T
    z, H = modo(P)
    cov_b = cho_solve(cho_factor(H), np.vstack([np.eye(p), np.zeros((n, p))]))[:p]
    return {"sigma2": sig ** 2, "phi": phi, "frontera": bool(phi > 0.98 or phi < 0.02), "beta": z[:p].tolist(),
            "cov_beta": cov_b.tolist(), "loglik_marg": float(-r.fun), "conv": bool(r.success), "theta": r.x.tolist(),
            "u": z[p:].tolist()}
