# Revisión de árbitro del artículo A3

Registro de la revisión por pares interna: qué se objetó, qué cambió, qué se mantuvo y por qué. Antes estaba como anexo B del manuscrito; se movió aquí para que el manuscrito que se envía a la revista no lo lleve. El código del recálculo del árbitro sirvió de base para `pruebas/test_calculos.py`.

## Contenido

1. Primera ronda: arbitraje simulado (29 de octubre de 2026)
2. Segunda ronda: reestructuración como repositorio de replicación (30 de octubre de 2026)
3. Tercera ronda: auditoría de determinismo (30 de octubre de 2026)

## 1. Primera ronda: arbitraje simulado (29 de octubre de 2026)
A otra instancia de Claude (Anthropic), la misma herramienta usada para preparar el artículo, se le pidió actuar como árbitro; no es una revisión por pares externa. Dictamen recibido el 29 de octubre de 2026 (revisión mayor). El árbitro recalculó las cifras principales y todas coincidieron; sus objeciones son de interpretación. La tesis central se debilitó: el título, el resumen y las conclusiones lo dicen.

| Objeción | Qué cambió | Qué se mantuvo y por qué |
|---|---|---|
| B1. La caída de la I no tiene línea base; con ruido espacial la I cae 10 % a 16 %; en hipertensivas la caída real no supera al placebo; la I es una razón. | Placebo calibrado a la I de Moran de cada covariable (200 réplicas), por bloque y en conjunto; exceso y valor p (tabla 2, figuras 2 y 3); numerador y denominador (sección 5.2); valor esperado bajo H0. Título, resumen, introducción y secciones 5.2, 7.1 y 8 reescritos. Se retiró que la certificación «determina buena porción del mapa hipertensivo». | La descomposición por bloques, ahora como exceso sobre el placebo. |
| M1. Shapley no monótono, bloques de tamaño distinto, sin incertidumbre. | Se reportan los 15 pasos no monótonos; el placebo tiene el tamaño y la estructura de cada bloque; intervalos de placebo. «16 % frente a 13 %» se reemplazó por los excesos (14,2 y 8,7) con sus intervalos. | |
| M2. Elegir contagio por LM y AIC es identificación muy débil. | Se quitaron «coherente con difusión», «forma de difusión» y la implicación sobre corredores; se añadió un panel anual con rezago espacio-temporal (tabla 3), que también muestra co-movimiento en el cáncer de próstata. | Las pruebas LM y el SDM/SDEM, como descripción de la forma de la dependencia. |
| M3. El suavizado global atenúa de forma diferencial; heteroscedasticidad; con pesos cambian los efectos directos. | Dos especificaciones coprincipales con pesos de precisión o conteos (SLX ponderado con errores agrupados; Poisson-BYM2 con rezagos); criterio de supervivencia; estimando único, el efecto total. | El SDEM sin pesos, como referencia. |
| M4. «Las perinatales no se agrupan entre municipios» es falso según el modelo de conteo. | Retirado; la sección 5.1 dice que la I baja refleja ruido y que el BYM2 atribuye 26 % y 50 % de la varianza a estructura espacial. | |
| M5. PDET: el total es «el municipio y todos sus vecinos»; selección endógena; errores agrupados; completitud sin variación suficiente. | Lenguaje no causal; errores agrupados y de Conley; el efecto propio de PDET no se distingue de cero; la completitud se reporta solo como efecto total y se declara su limitación. | La asociación total de PDET, que sobrevive en ambas especificaciones. |
| M6. Kodaman et al. no hablan de altitud; el Chow no tiene potencia. | La altitud se presenta como marcador geográfico; se retiró la frase sobre homogeneidad entre regiones. | |
| M7. I10 está en el numerador de los códigos basura. | Indicador recalculado sin I10–I15: la asociación con la mortalidad hipertensiva desaparece, y se dice. | |
| M8. Leroux en la frontera; φ no es fracción de varianza sin escalar; sin WX. | BYM2 escalado con rezagos, por Laplace; casos de frontera reportados; limitación declarada (sin INLA ni MCMC). | |
| m1. Selección por pretest. | Declarada en la introducción y en la sección 3.1. | |
| m2. BH separado para directos y totales. | Estimando único (efecto total) y BH conjunto por especificación sobre las 224 combinaciones. | |
| m3. «Cota inferior/superior» y E[I] = −0,025. | Retirado; se reporta el valor esperado bajo H0. | |
| m4. Figuras en grises; anonimato. | Figuras rediseñadas para grises y enlazadas desde `figuras/gris/`; pies sin colores; autocitas y menciones marcadas con ‡. | Versiones en color en `figuras/`. |
| m5. Efecto de excluir San Andrés. | Declarado en la sección 3.1. | |
| m6. Enfoque de Spijker et al. | Discutido en la sección 2.1. | |

## 2. Segunda ronda: reestructuración como repositorio de replicación (30 de octubre de 2026)

Al convertir el artículo en un repositorio autosuficiente, todas las etapas se volvieron a correr desde los agregados de `data/intermedios/`, y estos se reconstruyeron desde los microdatos con el paquete común de Estadísticas Vitales (`scripts/insumos_eevv/`). Hubo tres cambios.

| Qué | Cambio | Efecto en el manuscrito |
|---|---|---|
| Selección de las 12 causas (pretest) | Se recalculó la I de Moran departamental de 2022 de las 14 series con la armonización de edades del paquete común y la contigüidad departamental derivada del MGN 2018 (`scripts/11_departamental.py`). Malformaciones congénitas, otras afecciones perinatales y dificultad respiratoria del recién nacido dan −0,09, −0,11 y 0,08 (p de 0,69, 0,76 y 0,12), frente a 0,37, 0,31 y 0,27 en el trabajo previo. Con la misma matriz, las razones del trabajo previo dan 0,34, 0,31 y 0,25: la diferencia está en cómo ese trabajo agrupaba las edades de los menores de un año, no en la matriz. Las otras nueve causas conservan p < 0,05. | Secciones 1, 3.1, 5.1 y 7.2: se informa que tres causas de la selección no se agrupan entre departamentos con una agrupación de edades correcta. Los resultados municipales no cambian. El artículo descriptivo previo debe corregirse (pendiente fuera de este repositorio). |
| Leroux sin escalar de la versión 1 | Se retiró de la etapa de modelos (lo reemplazó el BYM2 escalado en la primera ronda). | Ninguno: el manuscrito ya no lo citaba. |
| Permutaciones de esda | Ahora usan una semilla fija (`np.random.seed(20261024)` antes de cada prueba), para que dos corridas den los mismos valores p. | La nota de la tabla 2 decía que la I inicial de las malformaciones tenía p = 0,002 con permutaciones sin semilla; con la semilla fija es 0,001, como las demás series, y la nota se corrigió. |

## 3. Tercera ronda: auditoría de determinismo (30 de octubre de 2026)

Una auditoría de replicación, hecha por otra instancia de Claude, clonó el primer commit y pidió las protecciones contra el no determinismo que aparecieron en A1, A2 y A4. Se aplicaron y el nivel A completo se corrió dos veces en limpio con semillas de hash distintas (`docs/REPLICACION.md`, sección 5).

| Qué | Cambio | Efecto en el manuscrito |
|---|---|---|
| Hilos y semilla de hash | El Makefile y `scripts/_entorno.py` fijan un hilo en BLAS, numba y polars y `PYTHONHASHSEED=0`; el script se relanza si faltan. | Ninguno. |
| Sumas de los esperados | Las sumas de reales por grupo (E municipal y departamental, esperados de la completitud) se hacen con numpy después de ordenar (`datos.suma_ordenada`), no con `group_by` de polars, cuyo orden de suma depende de los hilos. Los agregados versionados se regeneraron: cambian 1 o 2 unidades en el último dígito de E (diferencia relativa máxima 4e-16). | Ninguno: 202 comprobaciones coinciden. |
| Uniones | Las uniones de polars declaran `maintain_order="left"` y el panel se ordena por municipio y año antes de estimar. | Ninguno. |
| Panel de la leucemia | El suavizado anual contrae todas las SMR al promedio en los ocho años; la regresión quedaba sobre ruido de redondeo y daba γ = 1 con error estándar 1e-14, un número que cambiaba con el último dígito de E. `feols_dos_vias` ahora declara no estimable un ajuste sin variación dentro de los efectos fijos. En las otras afecciones perinatales la contracción ocurre en seis de los ocho años y δ es estimable pero prácticamente cero. | Nota de la tabla 3 precisada: la leucemia no es estimable y las afecciones perinatales dan δ menor que 0,0001. |

