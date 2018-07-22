"""go-basic.obo: fmt(1.2) rel(2018-07-16) 47,267 GO Terms; optional_attrs(relationship)"""

# Versions:
#    go-basic.obo: fmt(1.2) rel(2018-07-16) 47,267 GO Terms; optional_attrs(relationship)
#    goslim_generic.obo: fmt(1.2) rel(None) 232 GO Terms

# pylint: disable=line-too-long
SECTIONS = [ # 56 GO IDs placed into 13 sections; 0 unplaced GO IDs
    # ("New Section", [
    # ]),
    ("microtubule", [ # 3 GO-headers
        "GO:0007017",   # BP *    10 uGOs   298   8 L02 D02 R02 A     .... .r.. microtubule-based process
        "GO:0005856",   # CC **    2 uGOs   273   9 L04 D05 R07 ABD   .... p... cytoskeleton
        "GO:0005874",   # CC **    3 uGOs    29   7 L05 D06 R10 ABDEH P... p... microtubule
    ]),
    ("ubiquitin", [ # 4 GO-headers
        "GO:0006511",   # BP **    1 uGOs    40   3 L07 D08 R09 AC    .... prdu ubiquitin-dependent protein catabolic process
        "GO:0000151",   # CC *     2 uGOs    50  10 L03 D04 R05 ABC   .... .... ubiquitin ligase complex
        "GO:0031371",   # CC **    2 uGOs     4   4 L03 D04 R05 ABC   .... .... ubiquitin conjugating enzyme complex
        "GO:0019787",   # MF **    2 uGOs    32  11 L03 D03 R03 A     .... .... ubiquitin-like protein transferase activity
    ]),
    ("supramolecular", [ # 1 GO-headers
        "GO:0099080",   # CC **    4 uGOs    99   5 L01 D01 R01 H     .... .... supramolecular complex
    ]),
    ("protein", [ # 10 GO-headers
        "GO:0019538",   # BP **    5 uGOs  2201  22 L03 D04 R04 C     .... .rdu protein metabolic process
        "GO:0006464",   # BP **    5 uGOs  1318  62 L05 D06 R06 AC    .... .rdu cellular protein modification process
        "GO:0006508",   # BP **    5 uGOs   269   6 L04 D05 R05 C     .... .rdu proteolysis
        "GO:0006468",   # BP **    1 uGOs   261  20 L06 D07 R07 AC    .... .rdu protein phosphorylation
        "GO:1905368",   # CC *     3 uGOs    38   8 L03 D03 R03 C     .... .... peptidase complex
        "GO:0140096",   # MF **    2 uGOs   443  58 L02 D02 R02 A     .... .... catalytic activity, acting on a protein
        "GO:0008233",   # MF **    1 uGOs    75   4 L03 D03 R03 A     .... .... peptidase activity
        "GO:0008276",   # MF **    2 uGOs    33   7 L03 D05 R05 A     .... .... protein methyltransferase activity
        "GO:0061134",   # MF **    5 uGOs    15   3 L03 D03 R03 D     .... .... peptidase regulator activity
        "GO:0031386",   # MF **    1 uGOs     0   0 L01 D01 R01 E     .... .... protein tag
    ]),
    ("DNA_damage", [ # 2 GO-headers
        "GO:0006974",   # BP **    1 uGOs   210  10 L04 D04 R04 AF    .... prdu cellular response to DNA damage stimulus
        "GO:0006281",   # BP **    3 uGOs   135  15 L05 D07 R07 ACF   .... prdu DNA repair
    ]),
    ("phosph", [ # 3 GO-headers
        "GO:0016772",   # MF **    5 uGOs   604  10 L03 D03 R03 A     .... .... transferase activity, transferring phosphorus-containing groups
        "GO:0016301",   # MF **    2 uGOs   355 133 L04 D04 R04 A     .... .... kinase activity
        "GO:0016791",   # MF *     2 uGOs   179  70 L05 D05 R05 A     .... .... phosphatase activity
    ]),
    ("seeds", [ # 3 GO-headers
        "GO:0010344",   # BP **    1 uGOs     0   0 L02 D03 R09 DEHLM P... .... seed oilbody biogenesis
        "GO:0010169",   # CC **    1 uGOs     0   0 L03 D04 R07 ABC   .... .... thioglucosidase complex
        "GO:0010180",   # MF **    1 uGOs     0   0 L04 D04 R04 B     .... .... thioglucosidase binding
    ]),
    ("resistance", [ # 5 GO-headers
        "GO:0009404",   # BP **    8 uGOs   118  11 L03 D03 R03 AC    .... .... toxin metabolic process
        "GO:0098754",   # BP **    3 uGOs    56   7 L01 D01 R04 FQ    P... .... detoxification
        "GO:0009753",   # BP **    1 uGOs     7   3 L04 D05 R05 F     .... .... response to jasmonic acid
        "GO:0009635",   # BP **    1 uGOs     6   3 L03 D04 R04 F     .... .... response to herbicide
        "GO:0009864",   # BP **    1 uGOs     0   0 L03 D09 R09 ABFIJK P... .... induced systemic resistance, jasmonic acid mediated signaling pathway
    ]),
    ("cell junction", [ # 1 GO-headers
        "GO:0034330",   # BP **    7 uGOs    58   4 L03 D03 R03 AH    .... .... cell junction organization
    ]),
    ("cell cycle", [ # 2 GO-headers
        "GO:0022403",   # BP **    2 uGOs    46  12 L02 D02 R02 P     .... .... cell cycle phase
        "GO:0009524",   # CC **    1 uGOs     0   0 L04 D04 R07 AB    .... .... phragmoplast
    ]),
    ("metabolic", [ # 5 GO-headers
        "GO:0005975",   # BP *     5 uGOs   804   9 L03 D03 R03 C     .... .rdu carbohydrate metabolic process
        "GO:0006520",   # BP *     7 uGOs   535  17 L03 D06 R06 AC    .... .rdu cellular amino acid metabolic process
        "GO:0006790",   # BP **   10 uGOs   419  63 L03 D03 R03 AC    .... .rdu sulfur compound metabolic process
        "GO:0009617",   # BP **    1 uGOs    89   3 L03 D05 R05 FI    .... p... response to bacterium
        "GO:0009809",   # BP **    1 uGOs     7   3 L05 D06 R06 AC    .... .r.. lignin biosynthetic process
    ]),
    ("cellular_component", [ # 8 GO-headers
        "GO:0031090",   # CC **    2 uGOs   336   8 L02 D02 R03 DEF   P... p... organelle membrane
        "GO:0031224",   # CC **    3 uGOs   314   8 L02 D02 R03 FG    .... p... intrinsic component of membrane
        "GO:0031410",   # CC *     1 uGOs   283  25 L04 D05 R07 ABD   .... p... cytoplasmic vesicle
        "GO:0098805",   # CC **    1 uGOs   205  18 L02 D02 R02 F     .... p... whole membrane
        "GO:0005783",   # CC **    1 uGOs   107   5 L04 D05 R07 ABD   P... p... endoplasmic reticulum
        "GO:0005739",   # CC *     2 uGOs    96   4 L04 D05 R07 ABD   .... p... mitochondrion
        "GO:0005773",   # CC *     5 uGOs    85   6 L04 D05 R07 ABD   .... p... vacuole
        "GO:0005811",   # CC **    2 uGOs     3   2 L04 D05 R07 ABD   .... .... lipid droplet
    ]),
    ("Broad", [ # 9 GO-headers
        "GO:0009987",   # BP **    1 uGOs 18691  55 L01 D01 R01 A     .... .rdu cellular process
        "GO:0008152",   # BP **   12 uGOs  9786  21 L01 D01 R01 C     .... .rdu metabolic process
        "GO:0071840",   # BP *     2 uGOs  3585   2 L01 D01 R01 H     .... .... cellular component organization or biogenesis
        "GO:0051704",   # BP *     1 uGOs  2374  16 L01 D01 R01 I     .... .rdu multi-organism process
        "GO:0009056",   # BP **    5 uGOs  1839  10 L02 D02 R02 C     .... .rdu catabolic process
        "GO:0019748",   # BP **    2 uGOs   573  40 L02 D02 R02 C     .... .r.. secondary metabolic process
        "GO:0044848",   # BP **    2 uGOs    72   5 L01 D01 R01 P     .... .... biological phase
        "GO:0044464",   # CC **    1 uGOs  3295 113 L01 D01 R02 AB    P... .... cell part
        "GO:0016020",   # CC *     1 uGOs  1047  15 L01 D01 R01 F     .... p... membrane
    ]),
]
