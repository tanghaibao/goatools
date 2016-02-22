
run_goea:
	python test_goea_local.py
	python test_goea_rpt_bonferroni.py goea_bonferroni.log
	python test_goea_statsmodels.py
	python test_get_paths.py
	python test_genes_cell_cycle.py
	python test_mapslim.py
	python test_ncbi_entrez_annotations.py
	python test_optional_fields.py
	python test_write_hier.py

run_goea_scripts:
	cd ..; make goea_all

clean:
	rm -f goea_*.tsv goea_*.xlsx goea_*.log
	rm -f cell_cycle_gos.log cell_cycle_genes.txt
	rm -f genes_NCBI_hsa_All.pyc

clobber:
	@make --no-print-directory clean
	rm -f go-basic.obo
