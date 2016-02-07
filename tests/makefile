
run_goea:
	python test_goea_local.py
	python test_goea_rpt_bonferroni.py goea_bonferroni.log
	python test_goea_rpt_fdr_bh.py goea_fdr_bh.log
	python test_goea_statsmodels.py

run_goea_scripts:
	cd ..; make goea_all

clean:
	rm -f goea_*.tsv goea_*.xlsx goea_*.log

clobber:
	@make --no-print-directory clean
	rm -f go-basic.obo
