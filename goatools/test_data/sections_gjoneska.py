"""Sections variable used for grouping Gjoneska 2015 GO IDs."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


SECTIONS = [ # 18 sections
    ("immune", [ # 15 GO-headers
        "GO:0002376", # BP    564 L01 D01 M     immune system process
        "GO:0002682", # BP  1,183 L02 D02 AB    regulation of immune system process
        "GO:0030155", # BP    246 L02 D02 AB    regulation of cell adhesion
        "GO:0006955", # BP    100 L02 D02 GM    immune response
        "GO:0001817", # BP    476 L03 D03 AB    regulation of cytokine production
        "GO:0001775", # BP    162 L03 D03 CD    cell activation
        "GO:0001816", # BP    110 L03 D03 DK    cytokine production
        "GO:1903037", # BP    155 L04 D04 AB    regulation of leukocyte cell-cell adhesion
        "GO:0034097", # BP     59 L04 D04 G     response to cytokine
        "GO:0006954", # BP     25 L04 D04 G     inflammatory response
        "GO:0045087", # BP     25 L03 D04 GM    innate immune response
        "GO:0002521", # BP     72 L05 D05 CDF   leukocyte differentiation
        "GO:0007229", # BP      0 L05 D05 AB    integrin-mediated signaling pathway
        "GO:0050900", # BP     57 L02 D06 CDMN  leukocyte migration
        "GO:0042130", # BP      9 L07 D08 AB    negative regulation of T cell proliferation
        #"GO:0002252", # BP    138 L02 D02 L     immune effector process
    ]),
    ("viral/bacteria", [ # 4 GO-headers
        "GO:0016032", # BP    301 L03 D04 CJ    viral process
        "GO:0050792", # BP    119 L03 D04 AB    regulation of viral process
        "GO:0098542", # BP     37 L03 D05 GJ    defense response to other organism
        "GO:0009617", # BP     12 L03 D05 GJ    response to bacterium
    ]),
    ("neuro", [ # 25 GO-headers
        "GO:0099531", # BP     32 L01 D01 U     presynaptic process in chemical synaptic Xmission
        "GO:0042391", # BP    117 L03 D03 A     regulation of membrane potential
        "GO:0050877", # BP     96 L03 D03 K     neurological system process
        "GO:0050808", # BP     20 L03 D03 CDI   synapse organization
        "GO:0007272", # BP     13 L03 D03 CD    ensheathment of neurons
        "GO:0051960", # BP    236 L04 D04 AB    regulation of nervous system development
        "GO:0050804", # BP    120 L03 D04 AB    modulation of synaptic transmission
        "GO:0097485", # BP     34 L04 D04 CD    neuron projection guidance
        "GO:0031644", # BP     30 L04 D04 AB    regulation of neurological system process
        "GO:0031175", # BP     14 L04 D04 CDI   neuron projection development
        "GO:0035418", # BP     14 L04 D04 H     protein localization to synapse
        "GO:0007399", # BP      0 L04 D04 F     nervous system development
        "GO:0050767", # BP    192 L05 D05 AB    regulation of neurogenesis
        "GO:0030182", # BP     71 L05 D05 CDF   neuron differentiation
        "GO:0099536", # BP     40 L04 D05 CDR   synaptic signaling
        "GO:0048666", # BP     29 L04 D05 CDF   neuron development
        "GO:0010001", # BP     17 L05 D05 CDF   glial cell differentiation
        "GO:0051969", # BP      5 L03 D05 AB    regulation of transmission of nerve impulse
        "GO:0022008", # BP      3 L05 D05 CDF   neurogenesis
        "GO:0007158", # BP      0 L04 D05 DP    neuron cell-cell adhesion
        "GO:0014002", # BP      1 L05 D06 CDF   astrocyte development
        "GO:0048812", # BP     27 L05 D07 CDFI  neuron projection morphogenesis
        "GO:0048667", # BP      6 L06 D07 CDFI  cell morphogenesis involved in neuron differen.
        "GO:0072578", # BP      5 L05 D07 CDHI  neurotransmitter-gated ion channel clustering
        "GO:0007409", # BP     23 L06 D08 CDFI  axonogenesis
    ]),
    ("cell death", [ # 6 GO-headers
        "GO:0010941", # BP    316 L03 D03 AB    regulation of cell death
        "GO:0008219", # BP    104 L03 D03 CD    cell death
        "GO:0060548", # BP    103 L04 D04 AB    negative regulation of cell death
        "GO:0097190", # BP     22 L04 D04 AB    apoptotic signaling pathway
        "GO:0097527", # BP      0 L04 D04 AB    necroptotic signaling pathway
        "GO:0008637", # BP      7 L05 D05 CI    apoptotic mitochondrial changes
    ]),
    ("lipid", [ # 7 GO-headers
        "GO:0006629", # BP    623 L03 D03 DE    lipid metabolic process
        "GO:0019216", # BP    243 L04 D04 AB    regulation of lipid metabolic process
        "GO:0032368", # BP    130 L04 D04 AB    regulation of lipid transport
        "GO:0033993", # BP    112 L04 D04 G     response to lipid
        "GO:0006869", # BP     93 L04 D05 DH    lipid transport
        "GO:0055088", # BP     10 L05 D05 A     lipid homeostasis
        "GO:0042158", # BP      3 L05 D06 CE    lipoprotein biosynthetic process
    ]),
    ("adhesion", [ # 3 GO-headers
        "GO:0022610", # BP    194 L01 D01 P     biological adhesion
        "GO:0030155", # BP    246 L02 D02 AB    regulation of cell adhesion
        "GO:0007155", # BP    165 L02 D02 P     cell adhesion
    ]),
    ("cell cycle", [ # 9 GO-headers
        "GO:0022402", # BP    463 L02 D02 C     cell cycle process
        "GO:0022403", # BP     46 L02 D02 S     cell cycle phase
        "GO:0051726", # BP    411 L03 D03 AB    regulation of cell cycle
        "GO:0051301", # BP     54 L03 D03 CD    cell division
        "GO:0007049", # BP     12 L03 D03 CD    cell cycle
        "GO:0070192", # BP     17 L03 D05 CIL   chromosome organization in meiotic cell cycle
        "GO:0007051", # BP     19 L03 D06 CDI   spindle organization
        "GO:0007067", # BP      1 L04 D06 CI    mitotic nuclear division
        "GO:0030071", # BP     11 L06 D09 AB    regulation of mitotic metaphase/anaphase transition
    ]),
    ("chromosome", [ # 9 GO-headers
        "GO:0032259", # BP    119 L02 D02 E     methylation
        "GO:0051983", # BP    108 L03 D03 AB    regulation of chromosome segregation
        "GO:0007059", # BP     11 L03 D03 CD    chromosome segregation
        "GO:0006325", # BP    184 L04 D04 CI    chromatin organization
        "GO:0051276", # BP    107 L04 D04 CI    chromosome organization
        "GO:0032204", # BP     29 L03 D06 AB    regulation of telomere maintenance
        "GO:0034502", # BP     21 L06 D06 H     protein localization to chromosome
        "GO:0031497", # BP     11 L05 D06 CI    chromatin assembly
        "GO:0006334", # BP      3 L06 D07 CI    nucleosome assembly
    ]),
    ("development", [ # 10 GO-headers
        "GO:0032502", # BP  3,173 L01 D01 F     developmental process
        "GO:0022414", # BP    847 L01 D01 L     reproductive process
        "GO:0050793", # BP  1,881 L02 D02 AB    regulation of developmental process
        "GO:0048856", # BP  1,016 L02 D02 F     anatomical structure development
        "GO:0048646", # BP    331 L02 D02 F     anatomical structure formation in morphogenesis
        "GO:0007568", # BP     18 L03 D03 DF    aging
        "GO:0022604", # BP    129 L04 D04 AB    regulation of cell morphogenesis
        "GO:0000902", # BP     65 L04 D05 CDFI  cell morphogenesis
        "GO:0045765", # BP     14 L04 D05 AB    regulation of angiogenesis
    ]),
    ("extracellular matrix", [ # 1 GO-headers
        "GO:0030198", # BP     27 L04 D04 CDI   extracellular matrix organization
    ]),
    ("ion", [ # 3 GO-headers
        "GO:0006811", # BP    422 L04 D04 H     ion transport
        "GO:0055085", # BP    330 L04 D04 H     transmembrane transport
        "GO:0006874", # BP     33 L08 D09 ACD   cellular calcium ion homeostasis
    ]),
    ("localization", [ # 3 GO-headers
        "GO:0051179", # BP  2,142 L01 D01 H     localization
        "GO:0040011", # BP    394 L01 D01 N     locomotion
        "GO:0032879", # BP  1,682 L02 D02 AB    regulation of localization
    ]),
    ("membrane", [ # 1 GO-headers
        "GO:0061024", # BP    273 L03 D03 CI    membrane organization
    ]),
    ("metabolic", [ # 7 GO-headers
        "GO:0008152", # BP  6,418 L01 D01 E     metabolic process
        "GO:0019222", # BP  3,243 L02 D02 AB    regulation of metabolic process
        "GO:0009056", # BP  1,369 L02 D02 E     catabolic process
        "GO:0044281", # BP  2,139 L03 D03 DE    small molecule metabolic process
        "GO:0050790", # BP    620 L03 D03 A     regulation of catalytic activity
        "GO:0051186", # BP    373 L03 D03 CE    cofactor metabolic process
        "GO:0006259", # BP    300 L04 D06 CE    DNA metabolic process
    ]),
    ("phosphorylation", [ # 3 GO-headers
        "GO:0006793", # BP    798 L03 D03 CE    phosphorus metabolic process
        "GO:0016310", # BP    138 L05 D05 CE    phosphorylation
        "GO:0006468", # BP     97 L06 D07 CE    protein phosphorylation
    ]),
    ("signaling", [ # 4 GO-headers
        "GO:0023052", # BP    116 L01 D01 R     signaling
        "GO:0023051", # BP  1,364 L02 D02 AB    regulation of signaling
        "GO:0007165", # BP    717 L03 D03 AB    signal transduction
        "GO:0007267", # BP     99 L03 D04 CDR   cell-cell signaling
    ]),
    ("stimulus", [ # 4 GO-headers
        "GO:0050896", # BP  2,218 L01 D01 G     response to stimulus
        "GO:0048583", # BP  2,377 L02 D02 AB    regulation of response to stimulus
        "GO:0006950", # BP    492 L02 D02 G     response to stress
        "GO:0080134", # BP    940 L03 D03 AB    regulation of response to stress
    ]),
    ("prolif_differ", [ # 3 GO-headers
        "GO:0008283", # BP    158 L02 D02 D     cell proliferation
        "GO:0030154", # BP    494 L04 D04 CDF   cell differentiation
        "GO:0045595", # BP    828 L03 D03 AB    regulation of cell differentiation
        "GO:0042127", # BP    268 L03 D03 AB    regulation of cell proliferation
    ]),
]

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
