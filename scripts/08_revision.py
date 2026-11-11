#!/usr/bin/env python3
"""Nivel A. placebo, numerador y denominador, E[I], SLX ponderado y sin pesos con errores agrupados y de Conley, BH -> output/revision.json"""
import _entorno  # noqa: F401
from agrupamiento import modelos

if __name__ == "__main__":
    modelos.etapa_revision()
