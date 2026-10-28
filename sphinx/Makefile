# Build the GOATOOLS documentation site into _build/html

SPHINXBUILD  = sphinx-build
SPHINXAPIDOC = sphinx-apidoc
SPHINXOPTS   =

.PHONY: html apidoc clean

html: apidoc
	$(SPHINXBUILD) -b html $(SPHINXOPTS) . _build/html

# API pages generated from the package docstrings (test data excluded)
apidoc:
	$(SPHINXAPIDOC) -f -e -T -o _apidoc ../goatools ../goatools/test_data

clean:
	rm -rf _build _apidoc
