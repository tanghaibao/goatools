PYTHON := python


run:
	$(PYTHON) test_go_depth1
	$(PYTHON) test_get_unique_fields.py
	$(PYTHON) test_wr_py_goea_results.py
	$(PYTHON) test_altid_godag.py
	$(PYTHON) test_combine_nt_lists.py
	$(PYTHON) test_wr_tbl_subset.py
	$(PYTHON) test_pvalcalc.py
	$(PYTHON) test_go_print.py
	$(PYTHON) test_goea_local.py
	$(PYTHON) test_goea_rpt_bonferroni.py goea_bonferroni.log
	$(PYTHON) test_goea_statsmodels.py
	$(PYTHON) test_get_paths.py
	$(PYTHON) test_mapslim.py
	$(PYTHON) test_optional_fields.py
	$(PYTHON) test_write_hier.py
	$(PYTHON) test_annotations_gaf.py
	$(PYTHON) tests_methodsobj.py
	$(PYTHON) test_gaf_hdr.py
	$(PYTHON) test_go_draw.py
	$(PYTHON) test_read_gaf_allow_nd.py
	$(PYTHON) test_statsdescribe.py
	$(PYTHON) test_semantic_similarity.py
	$(PYTHON) similarity_top_parent.py noprt

dnld_dbs: go-basic.obo gene2go

# These tests require a download with wget
run_dnld:
	$(PYTHON) test_termcounts.py
	$(PYTHON) test_get_godag.py
	$(PYTHON) test_dnlds.py
	$(PYTHON) test_genes_cell_cycle.py
	$(PYTHON) test_ncbi_entrez_annotations.py
	$(PYTHON) test_write_summary_cnts.py
	$(PYTHON) test_nbt3102.py

run_goea_scripts:
	cd ..; make goea_all

go-basic.obo:
	wget http://geneontology.org/ontology/go-basic.obo

gene2go:
	wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
	gunzip gene2go.gz

dnld_obo:
	cd ..; make dnld_obo

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

clobber_dbs:
	rm -r go-basic.obo gene2go

clobber:
	@make --no-print-directory clean
	rm -f go-basic.obo
	rm -f gene2go
	rm -f gene2go*.gz
	rm -f *.*.s*d*
