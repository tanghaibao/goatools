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

# if the gene ontology files don't exist, download them
dnld_obo: $(GO_OBO_FILE) $(GOSLIM_OBO_FILE)

$(GO_OBO_FILE):
	@echo "downloading GO file: $(GO_OBO_FILE)"
	wget -O $(GO_OBO_FILE) $(GO_OBO_DOWNLOAD)

$(GOSLIM_OBO_FILE):
	@echo "downloading GOslim file: $(GOSLIM_OBO_FILE)"
	wget -O $(GOSLIM_OBO_FILE) $(GOSLIM_OBO_DOWNLOAD)

pylint:
	@git status -uno | perl -ne 'if (/(\S+.py)/) {printf "echo $$1\npylint -r no %s\n", $$1}' | tee tmp_pylint
	chmod 755 tmp_pylint
	tmp_pylint

clean:
	make clean_pyc
	rm -f goea*.xlsx goea.tsv GO_lineage.png
	cd tests; make --no-print-directory clean
	rm -f *.xlsx *.tsv *.log
	rm -f nbt3102_*

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
	@make --no-print-directory clean
	rm -f $(GO_OBO_FILE) $(GOSLIM_OBO_FILE)
	rm -f goa_human.gaf gene_association.*
	rm -f gene2go gene2go.gz
	rm -f *.gaf *.gaf.*
	rm -f goslim_*.obo goslim_*.owl
	rm -f goa_*.gpi goa_*.gpa

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

NOSETESTS := \
    tests/test_gosubdag_relationships.py \
    tests/test_gosubdag_mk.py \
    tests/test_go_depth1.py \
    tests/test_i92_relationship_parentchild.py \
    tests/test_i96_goea_ncbi.py \
    tests/test_gosearch_emptydict.py \
    tests/test_find_enrichment_overlap.py \
    tests/test_find_enrichment_run.py \
    tests/test_david_nts.py \
    tests/test_get_parents.py \
    tests/test_get_children.py \
    tests/test_optional_attributes.py \
    tests/test_genes_cell_cycle.py \
    tests/test_gpad_read.py \
    tests/test_semantic_similarity.py \
    tests/test_goea_errors.py \
    tests/test_ncbi_entrez_annotations.py \
    tests/test_wr_tbl_subset.py \
    tests/test_goea_local.py \
    tests/test_write_hier.py \
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
test:
	py.test tests/ --log-file=pytest.log

# This Representative test subset is automatically run for all push requests using Travis-CI.
# Running a subset of tests prevents Travis-CI from timeing out.
test_travis_subset:
	py.test $(NOSETESTS)

# Copyright (C) 2010-2017. Haibao Tang et al. All rights reserved.
