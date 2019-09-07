# For working with semantic similarity algorithm from:
# Yang, Haixuan et al.
# Bioinformatics (2012)
# Improving GO semantic similarity measures by exploring the ontology beneath the terms and modelling uncertainty

png:
	../../../scripts/go_plot.py --obo=fig1a.obo -o fig1a.png --id2gos=fig1a.anno
	../../../scripts/go_plot.py --obo=fig1b.obo -o fig1b.png --id2gos=fig1b.anno
	../../../scripts/go_plot.py --obo=fig2a.obo -o fig2a.png --id2gos=fig2a.anno

