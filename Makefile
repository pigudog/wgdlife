PY ?= python

.PHONY: test smoke fetch prepare
test:
	$(PY) -m pytest -q

smoke:
	$(PY) examples/run_wgd_example.py

fetch:
	$(PY) preprocess/fetch_msk_impact_50k.py

prepare:
	$(PY) preprocess/prepare_msk_impact_50k.py
