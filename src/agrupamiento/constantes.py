"""Constantes del artículo: periodo, causas, bloques de covariables, semillas."""
from __future__ import annotations

ANIOS = list(range(2018, 2026))
SIN_PANDEMIA = [2018, 2019, 2022, 2023, 2024, 2025]
SEMILLA = 20261024            # efectos del SDM y del SAR, permutaciones de la I de Moran
SEMILLA_PLACEBO = SEMILLA + 1
FUSION = {27493: 27615}       # Nuevo Belén de Bajirá -> Riosucio
EXCLUIR = {88001, 88564}      # San Andrés y Providencia: sin vecinos terrestres

CAUSAS = {  # código lista 6/67 -> (etiqueta corta, grupo)
    "todas": ("Todas las causas", "total"),
    604: ("Alzheimer, Parkinson y otras del sistema nervioso", "demas"),
    201: ("Cáncer de estómago", "neoplasias"),
    613: ("Malformaciones congénitas", "demas"),
    204: ("Cáncer de pulmón", "neoplasias"),
    406: ("Otras afecciones perinatales", "perinatales"),
    108: ("Neumonía e infecciones respiratorias agudas", "transmisibles"),
    210: ("Cáncer de próstata", "neoplasias"),
    212: ("Leucemia", "neoplasias"),
    302: ("Enfermedades hipertensivas", "circulatorio"),
    505: ("Ahogamiento", "externas"),
    404: ("Dificultad respiratoria del recién nacido", "perinatales"),
    103: ("Enfermedades por vectores y rabia", "transmisibles"),
    512: ("Homicidios", "externas"),
}
CLAVES = list(CAUSAS)

BLOQUES = {
    "S": ["nbi", "rural", "ldens", "indigena", "afro", "internet", "superior"],
    "G": ["altitud", "ldist_capital", "ldist_uci", "lcamas", "pdet"],
    "Q": ["sin_asistencia", "basura", "casa", "completitud"],
}
COVS = BLOQUES["S"] + BLOQUES["G"] + BLOQUES["Q"]
NOMBRE_COV = {
    "nbi": "NBI 2018 (%)", "rural": "Población rural (%)", "ldens": "Densidad (log hab/km²)",
    "indigena": "Población indígena (%)", "afro": "Población afrocolombiana (%)",
    "internet": "Viviendas con internet (%)", "superior": "Educación superior (%)",
    "altitud": "Altitud (m)", "ldist_capital": "Distancia a la capital (log km)",
    "ldist_uci": "Distancia a UCI (log km)", "lcamas": "Camas por 1.000 hab. (log)", "pdet": "Municipio PDET",
    "sin_asistencia": "Muertes por enfermedad sin asistencia (%)", "basura": "Códigos basura (%)",
    "casa": "Muertes en casa (%)", "solo_interrog": "Causa solo por interrogatorio (%)",
    "completitud": "Completitud departamental DANE (%)",
}
MACRO = {**{d: "Andina" for d in (5, 11, 15, 17, 25, 41, 54, 63, 66, 68, 73)},
         **{d: "Caribe" for d in (8, 13, 20, 23, 44, 47, 70, 88)},
         **{d: "Pacífica" for d in (19, 27, 52, 76)},
         **{d: "Orinoquía" for d in (50, 81, 85, 99)},
         **{d: "Amazonía" for d in (18, 86, 91, 94, 95, 97)}}


def clave(c) -> str:
    return str(c)
