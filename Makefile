PYTHON ?= python

.PHONY: setup pipeline dashboard clean

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py --server.port 8501 --server.address 0.0.0.0 --server.headless true

clean:
	rm -rf cell_counts.db outputs
