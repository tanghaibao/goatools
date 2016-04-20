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
	$(PYTHON) test_nbt3102.py
	$(PYTHON) test_write_summary_cnts.py
	$(PYTHON) tests_methodsobj.py

run_goea_scripts:
	cd ..; make goea_all

clean:
	rm -f nbt3102*.txt nbt3102*.xlsx
	rm -f goea_*.tsv goea_*.xlsx goea_*.log
	rm -f cell_cycle_*.log cell_cycle_*.txt
	rm -f *.pyc
	rm -f gene_association.*
	rm -f *.png
	rm -f tmp*

clobber:
	@make --no-print-directory clean
	rm -f go-basic.obo
	rm -f gene2go
	rm -f gene2go.gz
