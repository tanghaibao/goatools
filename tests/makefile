PYTHON := python

run:
	$(PYTHON) test_get_unique_fields.py
	$(PYTHON) test_wr_py_goea_results.py
	$(PYTHON) test_altid_godag.py
	$(PYTHON) test_combine_nt_lists.py
	$(PYTHON) test_wr_tbl_subset.py
	$(PYTHON) test_pvalcalc.py
	$(PYTHON) test_go_print.py
	$(PYTHON) test_dnlds.py
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
	$(PYTHON) test_write_summary_cnts.py
	$(PYTHON) tests_methodsobj.py
	$(PYTHON) test_nbt3102.py

run_goea_scripts:
	cd ..; make goea_all

clean_tmp:
	rm -f tmp*

clean:
	rm -f nbt3102*.txt nbt3102*.xlsx
	rm -f goea_*.tsv goea_*.xlsx goea_*.log
	rm -f cell_cycle_*.log cell_cycle_*.txt
	rm -f *.pyc
	rm -f gene_association.*
	rm -f *.png *.eps *.tif
	rm -f tmp*
	rm -f *.tsv
	rm -f nbt3102_goea.py nbt3102_goea_alt.py

clobber:
	@make --no-print-directory clean
	rm -f go-basic.obo
	rm -f gene2go
	rm -f gene2go*.gz
