#!/usr/bin/env python3
"""Nivel A. panel anual con rezago espacio-temporal -> output/panel.json"""
import _entorno  # noqa: F401
from agrupamiento import modelos

if __name__ == "__main__":
    modelos.etapa_panel()
