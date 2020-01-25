# configuration for Gene Ontology files
GO_OBO_FILE=go-basic.obo
GOSLIM_OBO_FILE=goslim_generic.obo

GO_OBO_DOWNLOAD=http://geneontology.org/ontology/go-basic.obo
GOSLIM_OBO_DOWNLOAD=http://www.geneontology.org/ontology/subsets/goslim_generic.obo
GOEA_FILES = data/study data/population data/association

# Example GO: translation factor activity RNA binding
GO = GO:0008135

# -------------------------------------------------------------------------------
# ---- Run GOEA -----------------------------------------------------------------
# -------------------------------------------------------------------------------
goea: $(GO_OBO_FILE)
	python scripts/find_enrichment.py --pval=0.05 --indent $(GOEA_FILES) --outfile results.txt

# -------------------------------------------------------------------------------
# ---- Compare 2+ GOEAS ---------------------------------------------------------
# -------------------------------------------------------------------------------
compare_gos:
	python scripts/compare_gos.py \
	data/compare_gos/tat_gos_simple1.tsv \
	data/compare_gos/tat_gos_simple2.tsv \
	-s data/compare_gos/sections.txt

compare_gos_wr:
	python scripts/compare_gos.py \
	data/compare_gos/tat_gos_simple1.tsv \
	data/compare_gos/tat_gos_simple2.tsv \
	-s data/compare_gos/sections.txt \
	--xlsx=tat_gos_simple.xlsx \
	-o tat_gos_simple.txt

# ---------------------------------------------------------------------------------------
# ---- Grouping -------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
USR := usr
DATE := 2018_0720
SEC    := $(USR)$(DATE)_goea
SECGO  := goea.tsv
SECIN  := $(SEC)_sections_in.txt
SECOUT := $(SEC)_sections.txt
SECTXT := $(SEC)_grouped_gos.txt
SECPY := $(SEC)_sections.py
GO_PLOT := scripts/go_plot.py

grpsec:
	echo $(SEC)
	wr_sections $(SECGO) -i $(SECIN) -o $(SECOUT) --txt=$(SECTXT) --py=$(SECPY)

grpvim:
	vim -p $(SECGO) $(SEC)_sections_in.txt $(SEC)_sections.txt $(SEC)_grouped_gos.txt

grpmisc:
	perl -ne 'if (/(Misc|$(TXT))/) {print}' grouped_gos.txt

grpre:
	perl -ne 'if (/(SECTION|$(TXT))/) {print}' $(SECTXT) | tee gos_$(TXT)

grpplt:
	wr_sections $(SECGO) -i $(SECIN) -o $(SECOUT) --txt=$(SECTXT) --py=$(SECPY)
	perl -ne 'if (/(SECTION|$(TXT))/) {print}' $(SECTXT) | tee gos_$(TXT)
	$(GO_PLOT) -s $(SECOUT) -i gos_$(TXT) -o aa_$(TXT).png


# -------------------------------------------------------------------------------
# ---- Sphinx-Generated Documentation -------------------------------------------
# -------------------------------------------------------------------------------

# GENERATE A TEMPORARY LOCAL WORKING COPY of Sphinx docs when developing documentation.
#
# User your browser to view the temporary local html files located at:
#   <LOCAL_GIT_ROOT>/sphinx/_build/html/index.html
#
# Html files generated using "mkdocs_practice" should not be committed or pushed.
# Use the "mkdocs_live" or "gh-pages" make target to generate html docs
# which will be saved in the "goatools" "GitHub Pages".
# "GitHub Pages" are public webpages hosted and published through the goatools repository.
.PHONY: mkdocs rmdocs gh-pages
mkdocs_practice:
	make -C sphinx/ apidoc html

# REMOVE THE TEMPORARY LOCAL WORKING COPY of Sphinx docs after a session of developing documentation.
rmdocs_practice:
	make clean_docgen


# Update Live on-line GOATOOLS Documentation
mkdocs_live:
	make gh-pages


# GENERATE SPHINX HTML DOCS DISPLAYED ON-LINE.
#
# Make target, "gh-pages" does the following:
#   1. Switches from the "master" branch to the "gh-pages" branch
#   2. While in the "gh-pages" branch, checks out from the "master" branch:
#      a. Sphinx control files
#      b. GOATOOLS source code, which contains docstrings
#         The docstrings are used to create Sphinx html documentation
#   3. Creates html documentation using Sphinx in the "gh-pages" branch.
#   4. Commits and pushes html docs from "gh-pages" branch.
#   5. Switches back to the master branch
gh-pages:
	git checkout gh-pages
	git rm -rf .
	git clean -dxf
	git checkout HEAD .nojekyll
	git checkout master sphinx goatools scripts
	make -C sphinx/ apidoc html
	mv -fv sphinx/_build/html/* .
	mv -fv _apidoc/* .
	rm -rf sphinx/ goatools/
	git add -A
	git commit -m "Generated gh-pages for `git log master -1 --pretty=short --abbrev-commit`"
	git push
	git checkout master


# -------------------------------------------------------------------------------
# ---- Run Scripts --------------------------------------------------------------
# -------------------------------------------------------------------------------
goea_scipy_pval: $(GO_OBO_FILE)
	python scripts/find_enrichment.py --pval=0.05 --indent $(GOEA_FILES) --pvalcalc fisher_scipy_stats

goea_basic: $(GO_OBO_FILE)
	python scripts/find_enrichment.py $(GOEA_FILES)

goea_xlsx: $(GO_OBO_FILE)
	python scripts/find_enrichment.py --pval=0.05 --indent $(GOEA_FILES) --outfile=goea.xlsx

goea_xlsx_bonferroni: $(GO_OBO_FILE)
	python scripts/find_enrichment.py --pval=0.05 --indent $(GOEA_FILES) --method=bonferroni --outfile=goea_bonferroni.xlsx

goea_tsv: $(GO_OBO_FILE)
	python scripts/find_enrichment.py --pval=0.05 --indent $(GOEA_FILES) --outfile=goea.tsv

goea_files: $(GO_OBO_FILE)
	python scripts/find_enrichment.py --pval=0.05 --indent $(GOEA_FILES) --outfile=goea.tsv,goea.xlsx

plot_go_pygraphviz: $(GO_OBO_FILE)
	python scripts/plot_go_term.py --term=$(GO)

plot_go_pydot: $(GO_OBO_FILE)
	python scripts/plot_go_term.py --term=$(GO) --engine=pydot

map_slim: $(GO_OBO_FILE) $(GOSLIM_OBO_FILE)
	python scripts/map_to_slim.py --association_file=data/association --slim_out=direct $(GO_OBO_FILE) $(GOSLIM_OBO_FILE)


goea_all: goea goea_basic goea_xlsx goea_xlsx_bonferroni goea_tsv goea_files

vim_compare_gos:
	vim -p \
	scripts/compare_gos.py \
	tests/test_compare_gos.py \
	tests/test_sorter_desc2nts.py \
	goatools/cli/compare_gos.py \
	goatools/cli/grouped.py \
	goatools/grouper/wrxlsx.py \
	goatools/grouper/sorter.py \
	goatools/grouper/sorter_gos.py \
	goatools/grouper/sorter_nts.py

# ./tests/test_optional_attributes.py
vim_attr:
	vim -p \
	./tests/test_optional_relationship.py \
	./goatools/obo_parser.py \
	./goatools/godag/obo_optional_attributes.py \
	./goatools/godag/typedef.py \
	./goatools/test_data/optional_attrs.py \
	./goatools/test_data/godag_timed.py

vim_ext:
	vim -p \
	tests/test_gpad_read.py \
	goatools/parsers/gpad_reader.py \
	goatools/anno/extensions/extensions.py \
	goatools/anno/extensions/extension.py

log_ver:
	git log -u goatools/__init__.py | tee tmp

# if the gene ontology files don't exist, download them
dnld_obo: $(GO_OBO_FILE) $(GOSLIM_OBO_FILE)

$(GO_OBO_FILE):
	@echo "downloading GO file: $(GO_OBO_FILE)"
	wget -O $(GO_OBO_FILE) $(GO_OBO_DOWNLOAD)

$(GOSLIM_OBO_FILE):
	@echo "downloading GOslim file: $(GOSLIM_OBO_FILE)"
	wget -O $(GOSLIM_OBO_FILE) $(GOSLIM_OBO_DOWNLOAD)

dnld_gene2go:
	@echo "Downloading gene2go"
	#wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
	gunzip gene2go.gz

pylint:
	@git status -uno | perl -ne 'if (/(\S+.py)/) {printf "echo $$1\npylint -r no %s\n", $$1}' | tee tmp_pylint
	chmod 755 tmp_pylint
	tmp_pylint

ls_tests:
	find tests -name \*.py -printf "%T@ %Tc %p\n" | sort -n

clean:
	make clean_pyc
	rm -f tmp_pylint
	rm -f goea*.xlsx goea.tsv GO_lineage.png
	cd tests; make --no-print-directory clean
	rm -f *.xlsx *.tsv *.log
	rm -f nbt3102_*
	rm -f data/gaf/goa_human_illegal.gaf.log
	rm -f gogrp_*.txt
	rm -f tests/data/gaf_missingsym.mgi.log
	rm -f tmp_test_sections.py
	rm -f ?_sec?_hdr*.txt
	rm -f compare_gos_*.txt
	rm -f transient_increase_hdrgos.txt
	rm -f tmp_test_sections_out.txt
	rm -f *.err

clean_pyc:
	find . -name \*.pyc | xargs rm -f
	rm -f py*.*.st*p

# Removes local files in master branch generated using Sphinx
clean_docgen:
	rm -rf ./_apidoc/
	rm -rf ./_modules/
	rm -rf ./_sources/
	rm -rf ./_static/
	rm -f genindex.html
	rm -f index.html
	rm -f objects.inv
	rm -f py-modindex.html
	rm -f search.html
	rm -f searchindex.js
	rm -rf sphinx/_apidoc/
	rm -rf sphinx/_build/
	rm -f fetch_associations.html
	rm -f find_enrichment.html
	rm -f goatools.html
	rm -f map_to_slim.html
	rm -f modules.html
	rm -f plot_go_term.html
	rm -f write_hierarchy.html

clobber:
	rm -f ./notebooks/go-basic.obo
	@make --no-print-directory clean
	rm -f $(GO_OBO_FILE) $(GOSLIM_OBO_FILE)
	rm -f goa_human.gaf gene_association.*
	rm -f gene2go gene2go.gz
	rm -f *.gaf *.gaf.*
	rm -f goslim_*.obo goslim_*.owl
	rm -f goa_*.gpi goa_*.gpa
	rm -f *.png
	rm -f gos_*
	rm -f cell_cycle_genes_*.txt
	rm -f *.gpa.gz
	rm -f gaf-eco-mapping-derived.txt
	rm -f *.gpad
	rm -f tmp*
	rm -f ids_stu_*.txt ids_pop_*.txt
	rm -f notebooks/*gpad*

ANNO = gaf
ANNO_HTTP = http://current.geneontology.org/annotations
dnld_anno:
	wget $(ANNO_HTTP)/aspgd.$(ANNO).gz; gunzip aspgd.$(ANNO).gz
	wget $(ANNO_HTTP)/cgd.$(ANNO).gz; gunzip cgd.$(ANNO).gz
	wget $(ANNO_HTTP)/dictybase.$(ANNO).gz; gunzip dictybase.$(ANNO).gz
	wget $(ANNO_HTTP)/ecocyc.$(ANNO).gz; gunzip ecocyc.$(ANNO).gz
	wget $(ANNO_HTTP)/fb.$(ANNO).gz; gunzip fb.$(ANNO).gz
	wget $(ANNO_HTTP)/genedb_lmajor.$(ANNO).gz; gunzip genedb_lmajor.$(ANNO).gz
	wget $(ANNO_HTTP)/genedb_pfalciparum.$(ANNO).gz; gunzip genedb_pfalciparum.$(ANNO).gz
	wget $(ANNO_HTTP)/genedb_tbrucei.$(ANNO).gz; gunzip genedb_tbrucei.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_chicken.$(ANNO).gz; gunzip goa_chicken.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_chicken_complex.$(ANNO).gz; gunzip goa_chicken_complex.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_chicken_isoform.$(ANNO).gz; gunzip goa_chicken_isoform.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_chicken_rna.$(ANNO).gz; gunzip goa_chicken_rna.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_cow.$(ANNO).gz; gunzip goa_cow.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_cow_complex.$(ANNO).gz; gunzip goa_cow_complex.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_cow_isoform.$(ANNO).gz; gunzip goa_cow_isoform.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_cow_rna.$(ANNO).gz; gunzip goa_cow_rna.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_dog.$(ANNO).gz; gunzip goa_dog.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_dog_complex.$(ANNO).gz; gunzip goa_dog_complex.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_dog_isoform.$(ANNO).gz; gunzip goa_dog_isoform.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_dog_rna.$(ANNO).gz; gunzip goa_dog_rna.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_human.$(ANNO).gz; gunzip goa_human.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_human_complex.$(ANNO).gz; gunzip goa_human_complex.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_human_isoform.$(ANNO).gz; gunzip goa_human_isoform.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_human_rna.$(ANNO).gz; gunzip goa_human_rna.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_pig.$(ANNO).gz; gunzip goa_pig.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_pig_complex.$(ANNO).gz; gunzip goa_pig_complex.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_pig_isoform.$(ANNO).gz; gunzip goa_pig_isoform.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_pig_rna.$(ANNO).gz; gunzip goa_pig_rna.$(ANNO).gz
	wget $(ANNO_HTTP)/gonuts.$(ANNO).gz; gunzip gonuts.$(ANNO).gz
	wget $(ANNO_HTTP)/gramene_oryza.$(ANNO).gz; gunzip gramene_oryza.$(ANNO).gz
	wget $(ANNO_HTTP)/jcvi.$(ANNO).gz; gunzip jcvi.$(ANNO).gz
	wget $(ANNO_HTTP)/mgi.$(ANNO).gz; gunzip mgi.$(ANNO).gz
	wget $(ANNO_HTTP)/pamgo_atumefaciens.$(ANNO).gz; gunzip pamgo_atumefaciens.$(ANNO).gz
	wget $(ANNO_HTTP)/pamgo_ddadantii.$(ANNO).gz; gunzip pamgo_ddadantii.$(ANNO).gz
	wget $(ANNO_HTTP)/pamgo_mgrisea.$(ANNO).gz; gunzip pamgo_mgrisea.$(ANNO).gz
	wget $(ANNO_HTTP)/pamgo_oomycetes.$(ANNO).gz; gunzip pamgo_oomycetes.$(ANNO).gz
	wget $(ANNO_HTTP)/pombase.$(ANNO).gz; gunzip pombase.$(ANNO).gz
	wget $(ANNO_HTTP)/pseudocap.$(ANNO).gz; gunzip pseudocap.$(ANNO).gz
	wget $(ANNO_HTTP)/reactome.$(ANNO).gz; gunzip reactome.$(ANNO).gz
	wget $(ANNO_HTTP)/rgd.$(ANNO).gz; gunzip rgd.$(ANNO).gz
	wget $(ANNO_HTTP)/sgd.$(ANNO).gz; gunzip sgd.$(ANNO).gz
	wget $(ANNO_HTTP)/sgn.$(ANNO).gz; gunzip sgn.$(ANNO).gz
	wget $(ANNO_HTTP)/tair.$(ANNO).gz; gunzip tair.$(ANNO).gz
	wget $(ANNO_HTTP)/wb.$(ANNO).gz; gunzip wb.$(ANNO).gz
	wget $(ANNO_HTTP)/zfin.$(ANNO).gz; gunzip zfin.$(ANNO).gz

dnld_anno_uniprot:
	wget $(ANNO_HTTP)/goa_uniprot_all.$(ANNO).gz; gunzip goa_uniprot_all.$(ANNO).gz
	wget $(ANNO_HTTP)/goa_uniprot_all_noiea.$(ANNO).gz; gunzip goa_uniprot_all_noiea.$(ANNO).gz

# Tests which run longer and have much functionality covered by other tests
#    tests/test_annotations_gaf.py \

# TBD: Add these to NOSETEST after edits:
#    tests/test_nbt3102.py \
#    tests/test_optional_fields.py \
#    tests/test_dnlds.py \
#    tests/test_get_godag.py \
#    tests/test_gpad_dnld.py \
#
#    tests/test_plot_get_parents.py \
#    tests/test_plot_objgoearesults.py \
#    tests/plt_i86obo.py \
#    tests/test_gosubdag_rcntobj.py \
#    tests/test_go_draw_basic.py \

#    tests/test_gosubdag_children.py \
#    tests/test_find_enrichment_overlap.py \
#    tests/test_find_enrichment_run.py \
#    tests/test_study_zero.py \
#    tests/test_gosubdag_mk.py \
#    tests/test_gpad_read.py \
#    tests/test_quickgo_xml.py \

NOSETESTS := \
		tests/test_parents_ancestors.py \
    tests/test_rpt_gene2go_evidencecodes.py \
    tests/test_sorter_sections.py \
    tests/test_sorter_desc2nts.py \
    tests/test_compare_gos.py \
    tests/test_wr_sections_txt.py \
    tests/test_altid_gosubdag.py \
    tests/test_dcnt_r01.py \
    tests/test_grprobj.py \
    tests/test_grpr_get_sections_2d.py \
    tests/test_sorter.py \
    tests/test_gosubdag_relationships.py \
    tests/test_go_depth1.py \
    tests/test_i92_relationship_parentchild.py \
    tests/test_i96_goea_ncbi.py \
    tests/test_gosearch_emptydict.py \
    tests/test_david_nts.py \
    tests/test_get_parents.py \
    tests/test_get_children.py \
    tests/test_optional_attributes.py \
    tests/test_genes_cell_cycle.py \
    tests/test_semantic_similarity.py \
    tests/test_goea_errors.py \
    tests/test_ncbi_entrez_annotations.py \
    tests/test_wr_tbl_subset.py \
    tests/test_goea_local.py \
    tests/test_write_hier.py \
		tests/test_cli_write_hierarchy.py \
    tests/test_go_print.py \
    tests/test_read_gaf_allow_nd.py \
    tests/test_write_summary_cnts.py \
    tests/test_pvalcalc.py \
    tests/test_altid_godag.py \
    tests/test_combine_nt_lists.py \
    tests/test_get_paths.py \
    tests/test_get_unique_fields.py \
    tests/test_go_draw.py \
    tests/test_goea_statsmodels.py \
    tests/test_goea_rpt_bonferroni.py \
    tests/test_wr_py_goea_results.py \
    tests/test_mapslim.py

# Run all tests. If you are submitting a pull request, all tests must pass.
pytest:
	#py.test -v tests/
	# make test_scripts
	python3 -m pytest -v tests | tee pytest.log
	# make chk_parsers
	# py.test tests/ --log-file=pytest.log
	grep FAIL pytest.log

test_scripts:
	python3 tests/test_cmds_goplot_rels.py all
	python3 tests/test_cmds_find_enrichment_md.py all
	python3 tests/test_cmds_plot_go.py

# Used to call GOATOOLS developers attention to illegal lines in parsed files
chk_parsers:
	tests/test_optional_attributes.py die

# This Representative test subset is automatically run for all push requests using Travis-CI.
# Running a subset of tests prevents Travis-CI from timeing out.
test_ci_subset:
	py.test --cov=goatools $(NOSETESTS)

# Copyright (C) 2010-2019. Haibao Tang et al. All rights reserved.
