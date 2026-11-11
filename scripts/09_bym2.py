#!/usr/bin/env python3
"""Nivel A. Poisson con efecto aleatorio BYM2 escalado y rezagos de las covariables, por Laplace -> output/bym2.json.
Sigue donde quedó si output/bym2.json existe; con --rehacer empieza de cero."""
import argparse

import _entorno  # noqa: F401
from agrupamiento import modelos

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rehacer", action="store_true")
    modelos.etapa_bym2("todo" if ap.parse_args().rehacer else None)
