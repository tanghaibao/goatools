PYTHON := python

run_goea:
	$(PYTHON) test_goea_local.py
	$(PYTHON) test_goea_rpt_bonferroni.py goea_bonferroni.log
	$(PYTHON) test_goea_statsmodels.py
	$(PYTHON) test_get_paths.py
	$(PYTHON) test_genes_cell_cycle.py
	$(PYTHON) test_mapslim.py
	$(PYTHON) test_optional_fields.py
	$(PYTHON) test_write_hier.py
	$(PYTHON) test_ncbi_entrez_annotations.py
	$(PYTHON) test_annotations_gaf.py

run_goea_scripts:
	cd ..; make goea_all

clean:
	rm -f goea_*.tsv goea_*.xlsx goea_*.log
	rm -f cell_cycle_gos.log cell_cycle_genes.txt
	rm -f genes_NCBI_hsa_All.pyc
	rm -f gene_association.*

clobber:
	@make --no-print-directory clean
	rm -f go-basic.obo
