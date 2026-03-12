# For working with semantic similarity algorithm from:
# Yang, Haixuan et al.
# Bioinformatics (2012)
# Improving GO semantic similarity measures by exploring the ontology beneath the terms and modelling uncertainty

png:
	python3 -m goatools go_plot --obo=fig1a.obo -o fig1a.png --id2gos=fig1a.anno
	python3 -m goatools go_plot --obo=fig1b.obo -o fig1b.png --id2gos=fig1b.anno
	python3 -m goatools go_plot --obo=fig2a.obo -o fig2a.png --id2gos=fig2a.anno
