#!/usr/bin/env python3
"""Nivel A. SDM y SDEM con efectos fijos de departamento -> output/resultados.json"""
import _entorno  # noqa: F401
from agrupamiento import modelos

if __name__ == "__main__":
    modelos.etapa_efectos_fijos()
