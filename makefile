# configuration for Gene Ontology files
GO_OBO_FILE=go-basic.obo
GOSLIM_OBO_FILE=goslim_generic.obo

GO_OBO_DOWNLOAD=http://geneontology.org/ontology/go-basic.obo
GOSLIM_OBO_DOWNLOAD=http://www.geneontology.org/ontology/subsets/goslim_generic.obo
GOEA_FILES = data/study data/population data/association

# Example GO: translation factor activity RNA binding
GO = GO:0008135


goea: $(GO_OBO_FILE)
	python scripts/find_enrichment.py --pval=0.05 --indent $(GOEA_FILES)

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


goea_all: goea goea_xlsx goea_xlsx_bonferroni goea_tsv goea_files


# if the gene ontology files don't exist, download them
$(GO_OBO_FILE):
	@echo "downloading GO file: $(GO_OBO_FILE)"
	wget -O $(GO_OBO_FILE) $(GO_OBO_DOWNLOAD)

$(GOSLIM_OBO_FILE):
	@echo "downloading GOslim file: $(GOSLIM_OBO_FILE)"
	wget -O $(GOSLIM_OBO_FILE) $(GOSLIM_OBO_DOWNLOAD)

clean:
	rm -f goea*.xlsx goea.tsv GO_lineage.png
	cd tests; make --no-print-directory clean

clobber:
	@make --no-print-directory clean
	rm -f $(GO_OBO_FILE) $(GOSLIM_OBO_FILE)
