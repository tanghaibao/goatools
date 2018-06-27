"""go-basic.obo: fmt(1.2) rel(2018-04-21) 47,216 GO Terms; optional_attrs(relationship)"""

# Versions:
#    go-basic.obo: fmt(1.2) rel(2018-04-21) 47,216 GO Terms; optional_attrs(relationship)
#    goslim_generic.obo: fmt(1.2) rel(None) 236 GO Terms

# pylint: disable=line-too-long
SECTIONS = [ # 122 GO IDs placed into 28 sections; 0 unplaced GO IDs
    # ("New Section", [
    # ]),
    ("immune", [ # 11 GO-headers
        "GO:0002376",   # BP **   59 uGOs  1796  19 L01 D01 R01 K     .... .rdu immune system process
        "GO:0001816",   # BP **   20 uGOs   678  52 L02 D02 R02 D     .... prdu cytokine production
        "GO:0045321",   # BP **   16 uGOs   443   6 L02 D03 R03 AK    .... .rdu leukocyte activation
        "GO:0045087",   # BP **    8 uGOs   288  13 L03 D04 R04 FK    .... prdu innate immune response
        "GO:0006954",   # BP **    4 uGOs   189   4 L04 D04 R04 F     .... prdu inflammatory response
        "GO:0050900",   # BP *     1 uGOs   157   9 L02 D05 R05 AGKN  .... .rdu leukocyte migration
        "GO:0002263",   # BP **    2 uGOs   149   4 L03 D03 R03 AFK   P... .... cell activation involved in immune response
        "GO:0007159",   # BP *     3 uGOs   125   3 L04 D04 R04 Q     .... .rdu leukocyte cell-cell adhesion
        "GO:0070661",   # BP *     6 uGOs    88   5 L02 D02 R02 O     .... .rdu leukocyte proliferation
        "GO:0071887",   # BP **    1 uGOs    75   5 L05 D05 R05 A     .... .rdu leukocyte apoptotic process
        "GO:0001909",   # BP *     5 uGOs    62   5 L02 D02 R02 S     .... .rdu leukocyte mediated cytotoxicity
    ]),
    ("viral/bacteria", [ # 5 GO-headers
        "GO:0016032",   # BP *     2 uGOs   361  47 L04 D04 R04 I     .... prdu viral process
        "GO:0098542",   # BP *     1 uGOs   149   8 L03 D05 R05 FI    .... p... defense response to other organism
        "GO:0009617",   # BP **    5 uGOs    89   3 L03 D05 R05 FI    .... p... response to bacterium
        "GO:0009620",   # BP **    1 uGOs    37   4 L03 D05 R05 FI    .... p... response to fungus
        "GO:0009615",   # BP **    3 uGOs    30   4 L03 D05 R05 FI    .... .... response to virus
    ]),
    ("neuro", [ # 16 GO-headers
        "GO:0007399",   # BP **    4 uGOs  1271   0 L04 D04 R05 DE    .... prdu nervous system development
        "GO:0022008",   # BP **    9 uGOs   705   2 L04 D04 R06 ADE   P... .rdu neurogenesis
        "GO:0007610",   # BP **   10 uGOs   280  26 L01 D01 R01 R     .... .rdu behavior
        "GO:0099536",   # BP *     6 uGOs   231   3 L03 D04 R04 AJ    .... p... synaptic signaling
        "GO:0031175",   # BP **    3 uGOs   223   7 L05 D05 R10 ADEH  P... prdu neuron projection development
        "GO:0050877",   # BP **   14 uGOs   217  14 L03 D03 R03 D     .... .rdu nervous system process
        "GO:0006836",   # BP **    3 uGOs   176  18 L04 D04 R04 G     .... .rdu neurotransmitter transport
        "GO:0048812",   # BP **    1 uGOs   155   5 L06 D07 R11 ADEH  P... p... neuron projection morphogenesis
        "GO:0048667",   # BP **    9 uGOs   149   4 L06 D06 R10 ADEH  P... p... cell morphogenesis involved in neuron differentiation
        "GO:0050808",   # BP **    8 uGOs   134   7 L03 D03 R03 AH    .... prd. synapse organization
        "GO:0061351",   # BP *     2 uGOs    75   8 L02 D02 R02 O     .... .rdu neural precursor cell proliferation
        "GO:0097479",   # BP *     3 uGOs    48   2 L05 D05 R05 G     .... p... synaptic vesicle localization
        "GO:0048167",   # BP **    7 uGOs    28  11 L03 D07 R09 ABJ   .... .... regulation of synaptic plasticity
        "GO:0035418",   # BP **    1 uGOs    22   5 L04 D04 R04 G     .... .r.u protein localization to synapse
        "GO:0019228",   # BP **    1 uGOs     3   0 L05 D05 R05 ABDJ  P... .rdu neuronal action potential
        "GO:0007158",   # BP **    1 uGOs     0   0 L04 D04 R04 Q     .... .... neuron cell-cell adhesion
    ]),
    ("cell death", [ # 4 GO-headers
        "GO:0008219",   # BP **   11 uGOs   503   4 L02 D02 R02 A     .... .rdu cell death
        "GO:0006915",   # BP **   12 uGOs   405  18 L04 D04 R04 A     .... prdu apoptotic process
        "GO:0001906",   # BP *     2 uGOs   113   3 L01 D01 R01 S     .... .rdu cell killing
        "GO:0008637",   # BP *     2 uGOs    16   4 L05 D05 R05 AH    P... .... apoptotic mitochondrial changes
    ]),
    ("lipid", [ # 8 GO-headers
        "GO:0006629",   # BP **   14 uGOs   948   8 L03 D03 R03 C     .... .rdu lipid metabolic process
        "GO:0010876",   # BP *     4 uGOs   260   1 L03 D03 R03 G     .... prdu lipid localization
        "GO:0033993",   # BP *     2 uGOs   196  35 L04 D04 R04 F     .... .... response to lipid
        "GO:0042157",   # BP *     1 uGOs    64   3 L04 D05 R05 C     .... .rdu lipoprotein metabolic process
        "GO:0097006",   # BP **    1 uGOs    47   0 L03 D03 R03 B     .... p... regulation of plasma lipoprotein particle levels
        "GO:0071825",   # BP *     2 uGOs    34   4 L04 D04 R04 AH    .... .... protein-lipid complex subunit organization
        "GO:0055088",   # BP **    1 uGOs    12   8 L05 D05 R05 B     .... .... lipid homeostasis
        "GO:0071402",   # BP *     1 uGOs    10   4 L03 D03 R03 AF    .... p... cellular response to lipoprotein particle stimulus
    ]),
    ("adhesion", [ # 3 GO-headers
        "GO:0022610",   # BP **    2 uGOs   350   4 L01 D01 R01 Q     .... .... biological adhesion
        "GO:0007155",   # BP **    6 uGOs   311   8 L02 D02 R02 Q     .... .rdu cell adhesion
        "GO:0031589",   # BP **    2 uGOs    44   7 L03 D03 R03 Q     .... .rdu cell-substrate adhesion
    ]),
    ("mitosis", [ # 4 GO-headers
        "GO:0000278",   # BP **   13 uGOs   363   5 L03 D03 R03 A     .... prdu mitotic cell cycle
        "GO:0140014",   # BP *    13 uGOs   123   1 L04 D06 R06 AH    .... prdu mitotic nuclear division
        "GO:0000281",   # BP **    1 uGOs    79   1 L04 D05 R06 A     .... prdu mitotic cytokinesis
        "GO:0098763",   # BP *     2 uGOs    12   3 L03 D03 R03 T     .... .... mitotic cell cycle phase
    ]),
    ("mieosis", [ # 3 GO-headers
        "GO:0051321",   # BP **    6 uGOs   301   1 L02 D03 R03 ALM   .... prdu meiotic cell cycle
        "GO:0045132",   # BP **    4 uGOs    90   5 L03 D04 R07 AHLM  P... p... meiotic chromosome segregation
        "GO:0051327",   # BP **    1 uGOs     0   0 L04 D04 R04 T     .... .... meiotic M phase
    ]),
    ("cell cycle", [ # 5 GO-headers
        "GO:0007049",   # BP **   14 uGOs   930   2 L02 D02 R02 A     .... prdu cell cycle
        "GO:0051301",   # BP **    3 uGOs   232  18 L02 D02 R02 A     .... prdu cell division
        "GO:0098813",   # BP *     7 uGOs   231   2 L03 D03 R04 A     .... p... nuclear chromosome segregation
        "GO:0044786",   # BP **    3 uGOs    80   3 L03 D09 R09 AC    .... p... cell cycle DNA replication
        "GO:0022403",   # BP **    3 uGOs    46  12 L02 D02 R02 T     .... .... cell cycle phase
    ]),
    ("damage/repair", [ # 2 GO-headers
        "GO:0006974",   # BP **   13 uGOs   210  10 L04 D04 R04 AF    .... prdu cellular response to DNA damage stimulus
        "GO:0006281",   # BP **    9 uGOs   135  15 L05 D07 R07 ACF   .... prdu DNA repair
    ]),
    ("chromatin", [ # 2 GO-headers
        "GO:0006325",   # BP **    9 uGOs   352  13 L03 D03 R05 AH    P... prdu chromatin organization
        "GO:0040029",   # BP **    4 uGOs   126   5 L06 D06 R06 BC    .... .... regulation of gene expression, epigenetic
    ]),
    ("DNA", [ # 6 GO-headers
        "GO:0051276",   # BP **    6 uGOs   682  18 L04 D04 R04 AH    .... prdu chromosome organization
        "GO:0097659",   # BP *     4 uGOs   614   2 L06 D08 R08 AC    .... .rdu nucleic acid-templated transcription
        "GO:0006259",   # BP **   23 uGOs   576  44 L04 D06 R06 AC    .... .rdu DNA metabolic process
        "GO:0007059",   # BP **    2 uGOs   239   1 L02 D02 R02 A     .... prdu chromosome segregation
        "GO:0006397",   # BP **    1 uGOs    66   7 L07 D08 R08 AC    .... prdu mRNA processing
        "GO:0051383",   # BP **    2 uGOs     9   1 L04 D05 R05 AH    .... .... kinetochore organization
    ]),
    ("development", [ # 6 GO-headers
        "GO:0032502",   # BP **    5 uGOs  6473  22 L01 D01 R01 E     .... .rdu developmental process
        "GO:0048856",   # BP **   26 uGOs  6065 195 L02 D02 R02 E     .... p... anatomical structure development
        "GO:0048646",   # BP *     1 uGOs   878 106 L02 D02 R04 E     P... .... anatomical structure formation involved in morphogenesis
        "GO:0008283",   # BP **    9 uGOs   492  31 L01 D01 R01 O     .... prdu cell proliferation
        "GO:0090066",   # BP *     1 uGOs   162   7 L03 D03 R03 B     .... .... regulation of anatomical structure size
        "GO:0048762",   # BP **    1 uGOs    76   8 L04 D04 R07 ADE   P... p... mesenchymal cell differentiation
    ]),
    ("extracellular matrix", [ # 1 GO-headers
        "GO:0030198",   # BP **    2 uGOs    51  10 L04 D04 R04 AH    .... prdu extracellular matrix organization
    ]),
    ("localization", [ # 7 GO-headers
        "GO:0051179",   # BP **   18 uGOs  4354  10 L01 D01 R01 G     .... .r.. localization
        "GO:0006810",   # BP **   12 uGOs  2996  34 L03 D03 R03 G     .... .rdu transport
        "GO:0006811",   # BP **   10 uGOs  1031   3 L04 D04 R04 G     .... .rdu ion transport
        "GO:0040011",   # BP **    2 uGOs   843   9 L01 D01 R01 N     .... .rdu locomotion
        "GO:0006928",   # BP **    8 uGOs   826   5 L02 D02 R02 A     .... .rdu movement of cell or subcellular component
        "GO:0032409",   # BP **    6 uGOs   134   5 L03 D03 R03 B     .... .... regulation of transporter activity
        "GO:0090130",   # BP **    2 uGOs    76   3 L02 D02 R02 D     .... .... tissue migration
    ]),
    ("membrane", [ # 4 GO-headers
        "GO:0055085",   # BP **    9 uGOs   740  26 L04 D04 R04 G     .... .rdu transmembrane transport
        "GO:0061024",   # BP *     2 uGOs   377  30 L03 D03 R03 AH    .... p... membrane organization
        "GO:0042391",   # BP **    4 uGOs   133  20 L03 D03 R03 B     .... .... regulation of membrane potential
        "GO:0033619",   # BP *     1 uGOs     9   3 L05 D06 R06 C     .... .... membrane protein proteolysis
    ]),
    ("metabolic", [ # 3 GO-headers
        "GO:0008152",   # BP **   16 uGOs  9805  21 L01 D01 R01 C     .... .rdu metabolic process
        "GO:0034641",   # BP **    7 uGOs  3598  23 L03 D03 R03 AC    .... .... cellular nitrogen compound metabolic process
        "GO:0006091",   # BP **    3 uGOs   217  10 L03 D03 R03 AC    .... .r.. generation of precursor metabolites and energy
    ]),
    ("phosphorylation", [ # 2 GO-headers
        "GO:0006793",   # BP **    7 uGOs  1372  28 L03 D03 R03 AC    .... .rdu phosphorus metabolic process
        "GO:0006464",   # BP **   13 uGOs  1319  62 L05 D06 R06 AC    .... .rdu cellular protein modification process
    ]),
    ("signaling", [ # 3 GO-headers
        "GO:0007154",   # BP **    3 uGOs  2452  14 L02 D02 R02 A     .... prdu cell communication
        "GO:0023052",   # BP **   30 uGOs  2310   7 L01 D01 R01 J     .... prdu signaling
        "GO:0007267",   # BP **    3 uGOs   656  31 L02 D03 R03 AJ    .... p... cell-cell signaling
    ]),
    ("stimulus", [ # 1 GO-headers
        "GO:0050896",   # BP **   31 uGOs  6004  16 L01 D01 R01 F     .... .rdu response to stimulus
    ]),
    ("vascular", [ # 5 GO-headers
        "GO:0072359",   # BP *     9 uGOs   640   0 L04 D04 R05 DE    .... p... circulatory system development
        "GO:0003013",   # BP **    4 uGOs   318   4 L03 D03 R03 D     .... .... circulatory system process
        "GO:0086001",   # BP *     1 uGOs    56   5 L05 D05 R05 B     .... pr.. cardiac muscle cell action potential
        "GO:0050817",   # BP **    1 uGOs    42   2 L02 D02 R02 D     .... .rdu coagulation
        "GO:0007596",   # BP **    5 uGOs    37   0 L03 D05 R05 BDF   P... prdu blood coagulation
    ]),
    ("cytoskeleton", [ # 5 GO-headers
        "GO:0006928",   # BP **    8 uGOs   826   5 L02 D02 R02 A     .... .rdu movement of cell or subcellular component
        "GO:0007010",   # BP **    1 uGOs   383  12 L04 D04 R04 AH    .... .rdu cytoskeleton organization
        "GO:0030030",   # BP **    3 uGOs   382   5 L03 D03 R03 AH    .... .rdu cell projection organization
        "GO:0007017",   # BP **   14 uGOs   297   8 L02 D02 R02 A     .... .r.. microtubule-based process
        "GO:0030029",   # BP **    4 uGOs   242   4 L02 D02 R02 A     .... .r.. actin filament-based process
    ]),
    ("reproduction", [ # 3 GO-headers
        "GO:0022414",   # BP **    2 uGOs  1219 106 L01 D01 R02 LM    P... .rdu reproductive process
        "GO:0009790",   # BP *     1 uGOs   649   2 L03 D04 R04 DE    .... prdu embryo development
        "GO:0009792",   # BP *     3 uGOs   170   2 L04 D05 R05 DE    .... p... embryo development ending in birth or egg hatching
    ]),
    ("vesicle", [ # 2 GO-headers
        "GO:0016192",   # BP **    9 uGOs   466  20 L04 D04 R04 G     .... pr.. vesicle-mediated transport
        "GO:0016050",   # BP **    1 uGOs   139  13 L04 D04 R04 AH    .... .... vesicle organization
    ]),
    ("catalytic", [ # 1 GO-headers
        "GO:0050790",   # BP *    14 uGOs   629  12 L03 D03 R03 B     .... .... regulation of catalytic activity
    ]),
    ("activation", [ # 1 GO-headers
        "GO:0001775",   # BP **    4 uGOs   496  13 L02 D02 R02 A     .... .rdu cell activation
    ]),
    ("hemostasis", [ # 5 GO-headers
        "GO:0042592",   # BP **   14 uGOs   474   6 L03 D03 R03 B     .... p... homeostatic process
        "GO:0050878",   # BP **    1 uGOs    97  10 L03 D03 R03 B     .... p... regulation of body fluid levels
        "GO:0042303",   # BP *     2 uGOs    56   3 L02 D02 R02 D     .... p... molting cycle
        "GO:0007599",   # BP **    1 uGOs    42   2 L04 D04 R04 B     .... .rdu hemostasis
        "GO:0001503",   # BP *     2 uGOs    32   6 L02 D02 R02 D     .... prdu ossification
    ]),
    ("broad", [ # 5 GO-headers
        "GO:0008150",   # BP **    1 uGOs 29625  29 L00 D00 R00       .... .rdu biological_process
        "GO:0009987",   # BP **    1 uGOs 18703  56 L01 D01 R01 A     .... .rdu cellular process
        "GO:0065007",   # BP **    8 uGOs 13064   3 L01 D01 R01 B     .... .... biological regulation
        "GO:0032501",   # BP **    6 uGOs  7544  68 L01 D01 R01 D     .... .rdu multicellular organismal process
        "GO:0009058",   # BP *     3 uGOs  3646  19 L02 D02 R02 C     .... .rdu biosynthetic process
    ]),
]
