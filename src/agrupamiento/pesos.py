"""Matrices de pesos espaciales.

- contigüidad tipo reina entre municipios, construida con los polígonos del MGN 2018 (DANE) generalizados a 0,0005 grados
  (`reina_desde_mgn`); se guarda como lista de adyacencia en data/intermedios/pesos_reina.csv y el nivel A la lee de ahí
  (`reina_desde_csv`), sin geopandas;
- 5 y 8 vecinos más cercanos e inversa de la distancia hasta 100 km, desde los puntos de referencia (`desde_puntos`);
- contigüidad entre departamentos, derivada de la municipal (`departamentos_desde_municipios`).
Todas se estandarizan por filas.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from . import rutas


def reina_desde_mgn(muns: list[int]) -> pl.DataFrame:
    """Lista de adyacencia (mun, vecino) por contigüidad tipo reina sobre el MGN 2018 (geometría fina)."""
    import geopandas as gpd
    import libpysal
    from shapely.validation import make_valid
    cods = [str(c).zfill(5) for c in muns]
    gdf = gpd.read_file(rutas.EXT63 / "mgn2018_mpio_fino.geojson")
    gdf["MPIO_CDPMP"] = gdf["MPIO_CDPMP"].astype(str).str.zfill(5)
    gdf = gdf.set_index("MPIO_CDPMP").loc[cods].reset_index()
    gdf["geometry"] = gdf.geometry.apply(make_valid)
    w = libpysal.weights.Queen.from_dataframe(gdf, ids=cods, use_index=False, silence_warnings=True)
    filas = [(int(i), int(j)) for i in cods for j in sorted(w.neighbors[i])]
    return pl.DataFrame(filas, schema={"mun": pl.Int32, "vecino": pl.Int32}, orient="row").sort(["mun", "vecino"])


def reina_desde_csv(muns: list[int]):
    import libpysal
    adj = pl.read_csv(rutas.INTER / "pesos_reina.csv")
    cods = [str(c).zfill(5) for c in muns]
    nb = {c: [] for c in cods}
    for i, j in adj.iter_rows():
        nb[str(i).zfill(5)].append(str(j).zfill(5))
    return libpysal.weights.W(nb, ids=cods, silence_warnings=True)


def desde_puntos(muns, lat, lon) -> dict:
    import libpysal
    import pyproj
    cods = [str(c).zfill(5) for c in muns]
    tr = pyproj.Transformer.from_crs(4326, 3116, always_xy=True)
    x, y = tr.transform(np.asarray(lon), np.asarray(lat))
    pts = np.column_stack([x, y]) / 1000.0   # km, MAGNA-SIRGAS / Colombia Bogotá
    W = {"knn5": libpysal.weights.KNN.from_array(pts, k=5, ids=cods),
         "knn8": libpysal.weights.KNN.from_array(pts, k=8, ids=cods)}
    umbral = 100.0
    D = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    nb, wt = {}, {}
    for i, c in enumerate(cods):
        js = np.where(D[i] <= umbral)[0]
        if len(js) == 0:      # sin vecinos a menos de 100 km: su vecino más cercano
            js = np.array([int(D[i].argmin())])
        nb[c] = [cods[j] for j in js]; wt[c] = list(1.0 / D[i, js])
    W["invdist"] = libpysal.weights.W(nb, wt, ids=cods, silence_warnings=True)
    forzados = int(sum(1 for i in range(len(cods)) if (D[i] <= umbral).sum() == 0))
    return W, {"umbral_km": umbral, "con_vecino_forzado": forzados}


def todas(cov: pl.DataFrame):
    """Las cuatro matrices, estandarizadas por filas, con su resumen."""
    muns = cov["mun"].to_list()
    W = {"reina": reina_desde_csv(muns)}
    Wp, extra = desde_puntos(muns, cov["lat"].to_numpy(), cov["lon"].to_numpy())
    W.update(Wp)
    info = {}
    for k, w in W.items():
        info[k] = {"islas": len(w.islands), "vecinos_media": float(w.mean_neighbors), "vecinos_min": int(w.min_neighbors),
                   "vecinos_max": int(w.max_neighbors)}
        w.transform = "r"
    info["invdist"].update(extra)
    return W, info


def simetrica_binaria(muns: list[int]) -> np.ndarray:
    """Matriz de contigüidad binaria y simétrica (para el efecto BYM2)."""
    adj = pl.read_csv(rutas.INTER / "pesos_reina.csv")
    pos = {m: i for i, m in enumerate(muns)}
    A = np.zeros((len(muns), len(muns)))
    for i, j in adj.iter_rows():
        A[pos[i], pos[j]] = 1.0
    return np.maximum(A, A.T)


def departamentos_desde_municipios(adj: pl.DataFrame) -> pl.DataFrame:
    """Dos departamentos son vecinos si algún municipio de uno toca a alguno del otro (reina municipal)."""
    adj = adj.with_columns((pl.col("mun") // 1000).alias("d1"), (pl.col("vecino") // 1000).alias("d2"))
    return adj.filter(pl.col("d1") != pl.col("d2")).select(pl.col("d1").alias("dp"), pl.col("d2").alias("vecino")).unique().sort(["dp", "vecino"])


def matriz_departamentos(dps: list[int]) -> np.ndarray:
    adj = pl.read_csv(rutas.INTER / "pesos_departamentos.csv")
    pos = {d: i for i, d in enumerate(dps)}
    W = np.zeros((len(dps), len(dps)))
    for i, j in adj.iter_rows():
        if i in pos and j in pos:
            W[pos[i], pos[j]] = 1.0
    return W / W.sum(1, keepdims=True)
