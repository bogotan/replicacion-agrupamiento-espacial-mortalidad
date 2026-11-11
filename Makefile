# Artículo A3. Todo corre con nice (máquina compartida). Ver README.md y docs/REPLICACION.md.
PY = nice -n 10 python3
# un hilo en BLAS, numba y polars: el número de hilos cambia el orden de las sumas y el último dígito de los reales
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1
export MKL_NUM_THREADS = 1
export VECLIB_MAXIMUM_THREADS = 1
export NUMEXPR_NUM_THREADS = 1
export NUMBA_NUM_THREADS = 1
export POLARS_MAX_THREADS = 1
# semilla de hash fija (scripts/_entorno.py; docs/REPLICACION.md, sección 5)
export PYTHONHASHSEED = 0
MEDIR = $(PY) scripts/insumos_eevv/medir.py tiempos_$@.tsv
NUEVOS ?= data/interim/63/intermedios_b

.PHONY: insumos verificar-insumos bases-eevv intermedios comparar nivel-b nivel-a nivel-a-rapido modelos pruebas pruebas-eevv todo limpiar

insumos:             ## baja lo que falte (externos del artículo y microdatos del DANE) y verifica huellas
	$(PY) scripts/00_descargar_insumos.py

verificar-insumos:   ## solo verifica huellas
	$(PY) scripts/00_descargar_insumos.py --verificar

bases-eevv:          ## paquete común: zips del DANE -> data/processed/ (unos 10 minutos)
	$(MEDIR) eevv_01 -- $(PY) scripts/insumos_eevv/01_extraer.py
	$(MEDIR) eevv_02 -- $(PY) scripts/insumos_eevv/02_base_defunciones.py
	$(MEDIR) eevv_03 -- $(PY) scripts/insumos_eevv/03_base_nacimientos.py
	$(MEDIR) eevv_04 -- $(PY) scripts/insumos_eevv/04_poblacion.py

intermedios:         ## etapas propias del nivel B: data/processed (o $$PROCESADO) + externos -> $(NUEVOS)
	$(MEDIR) 01_conteos -- $(PY) scripts/01_conteos.py
	$(MEDIR) 02_calidad -- $(PY) scripts/02_calidad_registro.py
	$(MEDIR) 03_pesos -- $(PY) scripts/03_pesos.py
	$(MEDIR) 04_intermedios -- $(PY) scripts/04_intermedios.py --salida $(NUEVOS)

comparar:            ## compara los intermedios reconstruidos con los versionados
	$(PY) scripts/05_comparar_intermedios.py --nuevo $(NUEVOS) --referencia data/intermedios

nivel-b: insumos bases-eevv pruebas-eevv intermedios comparar

modelos:             ## estimaciones desde data/intermedios/ (unos 22 minutos; el BYM2 es la parte lenta)
	$(MEDIR) 06_modelos -- $(PY) scripts/06_modelos.py
	$(MEDIR) 07_efectos_fijos -- $(PY) scripts/07_efectos_fijos.py
	$(MEDIR) 08_revision -- $(PY) scripts/08_revision.py
	$(MEDIR) 09_bym2 -- $(PY) scripts/09_bym2.py --rehacer
	$(MEDIR) 10_panel -- $(PY) scripts/10_panel.py
	$(MEDIR) 11_departamental -- $(PY) scripts/11_departamental.py

nivel-a: modelos     ## completo: estimaciones, figuras, tablas y verificación de cifras
	$(MEDIR) 12_figuras -- $(PY) scripts/12_figuras_tablas.py
	$(PY) pruebas/verificar_cifras.py

nivel-a-rapido:      ## usa output/*.json versionados (sin reestimar): figuras, tablas y verificación de cifras
	$(MEDIR) 12_figuras -- $(PY) scripts/12_figuras_tablas.py
	$(PY) pruebas/verificar_cifras.py

pruebas:             ## pruebas del artículo (insumos versionados, estimadores, cálculos) y cifras del texto
	$(PY) -m pytest pruebas/ -q -p no:cacheprovider --ignore=pruebas/insumos_eevv
	$(PY) pruebas/verificar_cifras.py

pruebas-eevv:        ## pruebas del paquete común (huellas de los zips, columnas, totales del DANE): requieren el nivel B
	$(PY) -m pytest pruebas/insumos_eevv -q -p no:cacheprovider


todo: nivel-b nivel-a pruebas

limpiar:
	rm -rf data/interim data/processed
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
