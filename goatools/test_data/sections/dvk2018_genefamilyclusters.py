"""go-basic.obo: fmt(1.2) rel(2018-05-07) 47,222 GO Terms; optional_attrs(relationship)"""

# Versions:
#    go-basic.obo: fmt(1.2) rel(2018-05-07) 47,222 GO Terms; optional_attrs(relationship)
#    goslim_generic.obo: fmt(1.2) rel(None) 235 GO Terms

# pylint: disable=line-too-long
SECTIONS = [ # 209 GO IDs placed into 71 sections; 0 unplaced GO IDs
    # ("New Section", [
    # ]),
    ("Immune", [ # 2 GO-headers
        "GO:0002376",   # BP *    23 uGOs  1796  19 L01 D01 R01 K     .... .rdu immune system process
        "GO:0006952",   # BP **   11 uGOs   663  11 L03 D03 R03 F     .... prdu defense response
    ]),
    ("Immune_adaptive", [ # 2 GO-headers
        "GO:0042110",   # BP *     5 uGOs   230   6 L04 D05 R05 AK    .... prdu T cell activation
        "GO:0050852",   # BP **    1 uGOs     4   0 L06 D10 R10 ABFJK .... prdu T cell receptor signaling pathway
    ]),
    ("Immune_antigen_presentation", [ # 9 GO-headers
        "GO:0019882",   # BP **   24 uGOs   128  14 L02 D02 R02 K     .... .rdu antigen processing and presentation
        "GO:0042611",   # CC **    2 uGOs     3   3 L04 D04 R07 ABCFG .... .... MHC protein complex
        "GO:0042824",   # CC *     1 uGOs     2   0 L03 D05 R13 ABCDEFG P... p... MHC class I peptide loading complex
        "GO:0042287",   # MF *     1 uGOs     9   3 L04 D04 R04 B     .... .... MHC protein binding
        "GO:0023023",   # MF **    2 uGOs     6   3 L03 D03 R03 B     .... .... MHC protein complex binding
        "GO:0046977",   # MF **    2 uGOs     2   2 L03 D03 R03 B     .... .... TAP binding
        "GO:0042605",   # MF **    1 uGOs     0   0 L03 D04 R04 B     .... .... peptide antigen binding
        "GO:0032395",   # MF **    1 uGOs     0   0 L04 D04 R04 D     .... .... MHC class II receptor activity
        "GO:0015433",   # MF **    1 uGOs     0   0 L06 D13 R13 AC    .... .... peptide antigen-transporting ATPase activity
    ]),
    ("Immune_innate", [ # 8 GO-headers
        "GO:0060429",   # BP *     2 uGOs  1099  20 L04 D04 R04 E     .... p... epithelium development
        "GO:0045087",   # BP *     5 uGOs   288  13 L03 D04 R04 FK    .... prdu innate immune response
        "GO:0008544",   # BP **    5 uGOs    93   1 L04 D04 R04 E     .... prdu epidermis development
        "GO:0050764",   # BP **    4 uGOs    24   7 L06 D07 R07 ABGH  .R.. .... regulation of phagocytosis
        "GO:0002468",   # BP *     2 uGOs    12   2 L03 D03 R03 K     .... prdu dendritic cell antigen processing and presentation
        "GO:0097048",   # BP *     2 uGOs     3   0 L06 D06 R06 A     .... .rdu dendritic cell apoptotic process
        "GO:0097026",   # BP *     1 uGOs     3   0 L06 D06 R06 AH    .... .rdu dendritic cell dendrite assembly
        "GO:0030280",   # MF **    1 uGOs     1   1 L02 D02 R02 G     .... .... structural constituent of epidermis
    ]),
    ("Immune_complement", [ # 3 GO-headers
        "GO:0006956",   # BP **    3 uGOs    29   4 L03 D07 R07 BCFK  .... .rdu complement activation
        "GO:0097278",   # BP *     2 uGOs     3   0 L02 D02 R02 Q     .... .rdu complement-dependent cytotoxicity
        "GO:0001848",   # MF **    1 uGOs     8   8 L03 D03 R03 B     .... .... complement binding
    ]),
    ("Immune_pattern_recog_rx", [ # 3 GO-headers
        "GO:0002221",   # BP **    5 uGOs   110   4 L05 D10 R10 ABFJK .... .... pattern recognition receptor signaling pathway
        "GO:0035325",   # MF **    2 uGOs     2   2 L04 D04 R04 B     .... .... Toll-like receptor binding
        "GO:0050786",   # MF **    1 uGOs     0   0 L04 D04 R04 B     .... .... RAGE receptor binding
    ]),
    ("differentiation_proliferation", [ # 2 GO-headers
        "GO:0045595",   # BP **    7 uGOs   921  41 L04 D04 R04 ABE   .R.. .... regulation of cell differentiation
        "GO:0008283",   # BP *     1 uGOs   492  31 L01 D01 R01 N     .... prdu cell proliferation
    ]),
    ("Bacteria/virus", [ # 7 GO-headers
        "GO:0044403",   # BP **    8 uGOs  1372   5 L03 D03 R03 I     .... pr.. symbiont process
        "GO:0035821",   # BP **    2 uGOs   789  15 L03 D03 R03 I     .... .... modification of morphology or physiology of other organism
        "GO:0051707",   # BP **    7 uGOs   547  14 L02 D04 R04 FI    .... .... response to other organism
        "GO:0009617",   # BP **    6 uGOs    89   3 L03 D05 R05 FI    .... p... response to bacterium
        "GO:1900561",   # BP **    2 uGOs     5   2 L03 D03 R03 C     .... .... dehydroaustinol metabolic process
        "GO:1900558",   # BP **    2 uGOs     5   2 L03 D03 R03 C     .... .... austinol metabolic process
        "GO:0008745",   # MF **    1 uGOs     0   0 L03 D05 R05 A     .... .... N-acetylmuramoyl-L-alanine amidase activity
    ]),
    ("cell_activation", [ # 1 GO-headers
        "GO:0050865",   # BP **    8 uGOs   273   8 L04 D04 R04 AB    .R.. .... regulation of cell activation
    ]),
    ("clotting", [ # 2 GO-headers
        "GO:0050817",   # BP **    2 uGOs    42   2 L02 D02 R02 D     .... .rdu coagulation
        "GO:0030219",   # BP *     2 uGOs     4   0 L05 D05 R09 ADEK  .... prdu megakaryocyte differentiation
    ]),
    ("detoxification", [ # 13 GO-headers
        "GO:0009410",   # BP *     2 uGOs   189  24 L03 D03 R03 F     .... .... response to xenobiotic stimulus
        "GO:0009698",   # BP **    2 uGOs    93  13 L03 D04 R04 AC    .... .r.. phenylpropanoid metabolic process
        "GO:0098754",   # BP **    6 uGOs    56   7 L01 D01 R04 FR    P... .... detoxification
        "GO:0009812",   # BP **    2 uGOs    47   9 L03 D03 R03 C     .... p... flavonoid metabolic process
        "GO:0046184",   # BP **    2 uGOs    24  12 L04 D04 R04 AC    .... .... aldehyde biosynthetic process
        "GO:0006063",   # BP **    1 uGOs    22   2 L04 D07 R07 AC    .... .... uronic acid metabolic process
        "GO:0019585",   # BP **    4 uGOs    15   5 L05 D08 R08 AC    .... .... glucuronate metabolic process
        "GO:1901685",   # BP **    2 uGOs     8   3 L04 D04 R04 AC    .... .... glutathione derivative metabolic process
        "GO:0016491",   # MF *     3 uGOs  2367  87 L02 D02 R02 A     .... .... oxidoreductase activity
        "GO:0015020",   # MF **    1 uGOs    12  12 L05 D05 R05 A     .... .... glucuronosyltransferase activity
        "GO:0004364",   # MF **    1 uGOs     0   0 L04 D04 R04 A     .... .... glutathione transferase activity
        "GO:0033906",   # MF **    1 uGOs     0   0 L05 D05 R05 A     .... .... hyaluronoglucuronidase activity
        "GO:0043295",   # MF **    1 uGOs     0   0 L03 D05 R05 B     .... .... glutathione binding
    ]),
    ("receptor", [ # 2 GO-headers
        "GO:0004888",   # MF **    2 uGOs   332  36 L03 D03 R03 D     .... .... transmembrane signaling receptor activity
        "GO:0030545",   # MF **    5 uGOs    46   4 L02 D02 R02 E     .... .... receptor regulator activity
    ]),
    ("Cytokine", [ # 6 GO-headers
        "GO:0001816",   # BP *     6 uGOs   678  52 L02 D02 R02 D     .... prdu cytokine production
        "GO:0034097",   # BP *    10 uGOs   185  30 L04 D04 R04 F     .... .rdu response to cytokine
        "GO:0000270",   # BP **    2 uGOs     8   3 L06 D06 R06 C     .... .... peptidoglycan metabolic process
        "GO:0005126",   # MF **    2 uGOs    92  48 L04 D04 R04 B     .... .... cytokine receptor binding
        "GO:0005125",   # MF **    2 uGOs     1   1 L04 D05 R05 BE    .... .... cytokine activity
        "GO:0016019",   # MF **    1 uGOs     0   0 L05 D05 R05 D     .... .... peptidoglycan receptor activity
    ]),
    ("vesicle", [ # 2 GO-headers
        "GO:0031982",   # CC **    2 uGOs   299   2 L03 D03 R03 D     .... p... vesicle
        "GO:0031410",   # CC *    11 uGOs   283  25 L04 D05 R07 ABD   .... p... cytoplasmic vesicle
    ]),
    ("Golgi", [ # 3 GO-headers
        "GO:1990668",   # BP **    1 uGOs     2   2 L06 D07 R07 AGH   .... .... vesicle fusion with endoplasmic reticulum-Golgi intermediate compartment (ERGIC) membrane
        "GO:0005794",   # CC *     3 uGOs    44   0 L04 D05 R07 ABD   P... p... Golgi apparatus
        "GO:0005798",   # CC *     2 uGOs    20   3 L05 D06 R08 ABD   .... p... Golgi-associated vesicle
    ]),
    ("endoplasmic_reticulum", [ # 2 GO-headers
        "GO:0046967",   # BP **    1 uGOs     0   0 L04 D05 R05 G     .... .... cytosol to ER transport
        "GO:0005783",   # CC **    6 uGOs   107   5 L04 D05 R07 ABD   P... p... endoplasmic reticulum
    ]),
    ("nucleus", [ # 4 GO-headers
        "GO:0005634",   # CC **    1 uGOs   510  12 L04 D05 R07 ABD   .... p... nucleus
        "GO:0044428",   # CC **    1 uGOs   489 125 L03 D04 R08 ABDE  P... .... nuclear part
        "GO:0005654",   # CC **    3 uGOs   126   0 L04 D05 R10 ABDEH P... p... nucleoplasm
        "GO:0048471",   # CC **    1 uGOs     7   1 L04 D04 R07 AB    .... p... perinuclear region of cytoplasm
    ]),
    ("extracellular", [ # 4 GO-headers
        "GO:0005576",   # CC **    1 uGOs   198   4 L01 D01 R01 I     .... p... extracellular region
        "GO:0044421",   # CC **    4 uGOs   193  39 L01 D01 R02 IJ    P... .... extracellular region part
        "GO:0005615",   # CC **    3 uGOs    50   0 L02 D02 R03 IJ    .... p... extracellular space
        "GO:0005198",   # MF **    1 uGOs    46  21 L01 D01 R01 G     .... .... structural molecule activity
    ]),
    ("ion", [ # 5 GO-headers
        "GO:0050801",   # BP **   13 uGOs   205   6 L05 D05 R05 B     .... .... ion homeostasis
        "GO:0010038",   # BP **   10 uGOs    82  19 L04 D04 R04 F     .... .... response to metal ion
        "GO:0051238",   # BP **    2 uGOs    12   5 L03 D04 R04 BG    .... .... sequestering of metal ion
        "GO:0043167",   # MF **    7 uGOs   206   2 L02 D02 R02 B     .... .... ion binding
        "GO:0031406",   # MF **    6 uGOs    61  17 L04 D04 R04 B     .... .... carboxylic acid binding
    ]),
    ("photo", [ # 5 GO-headers
        "GO:0009416",   # BP *     1 uGOs   117  13 L04 D04 R04 F     .... .... response to light stimulus
        "GO:0071482",   # BP **    2 uGOs    37   5 L05 D06 R06 AF    .... .... cellular response to light stimulus
        "GO:0009583",   # BP **    2 uGOs    20   4 L04 D05 R05 F     .... p... detection of light stimulus
        "GO:0097733",   # CC *     1 uGOs     4   0 L04 D07 R08 ABD   .... p... photoreceptor cell cilium
        "GO:0009881",   # MF **    2 uGOs     5   3 L03 D03 R03 D     .... .... photoreceptor activity
    ]),
    ("olfactory", [ # 1 GO-headers
        "GO:0004984",   # MF **    1 uGOs     1   1 L04 D04 R04 D     .... .... olfactory receptor activity
    ]),
    ("transcription", [ # 6 GO-headers
        "GO:0010468",   # BP **    6 uGOs   890  19 L05 D05 R05 BC    .R.. p... regulation of gene expression
        "GO:0097659",   # BP **   14 uGOs   617   2 L06 D08 R08 AC    .... .rdu nucleic acid-templated transcription
        "GO:0008023",   # CC **    2 uGOs    10   8 L02 D06 R12 ABCDEH .... .... transcription elongation factor complex
        "GO:0140110",   # MF **    1 uGOs    94   6 L01 D01 R01 F     .... .... transcription regulator activity
        "GO:0003700",   # MF **    4 uGOs    58  10 L02 D02 R02 F     .... .... DNA binding transcription factor activity
        "GO:0044212",   # MF **    6 uGOs    52  12 L05 D05 R05 B     .... .... transcription regulatory region DNA binding
    ]),
    ("Lipid", [ # 2 GO-headers
        "GO:0006629",   # BP *     6 uGOs   948   8 L03 D03 R03 C     .... .rdu lipid metabolic process
        "GO:0008289",   # MF **    3 uGOs    96  12 L02 D02 R02 B     .... .... lipid binding
    ]),
    ("Cell Death", [ # 3 GO-headers
        "GO:0008219",   # BP **    5 uGOs   503   4 L02 D02 R02 A     .... .rdu cell death
        "GO:0001906",   # BP **    4 uGOs   113   3 L01 D01 R01 Q     .... .rdu cell killing
        "GO:0005035",   # MF **    2 uGOs     3   2 L04 D04 R04 D     .... .... death receptor activity
    ]),
    ("Signaling", [ # 2 GO-headers
        "GO:0023052",   # BP *     6 uGOs  2307   7 L01 D01 R01 J     .... prdu signaling
        "GO:0007267",   # BP **    5 uGOs   657  31 L02 D03 R03 AJ    .... p... cell-cell signaling
    ]),
    ("chromatin", [ # 4 GO-headers
        "GO:0006325",   # BP **   11 uGOs   352  13 L03 D03 R05 AH    P... prdu chromatin organization
        "GO:0000785",   # CC **    4 uGOs    77   7 L04 D05 R09 ABDE  .... p... chromatin
        "GO:0003682",   # MF **    4 uGOs     8   4 L02 D02 R02 B     .... .... chromatin binding
        "GO:0042393",   # MF **    1 uGOs     6   4 L03 D03 R03 B     .... .... histone binding
    ]),
    ("protein_DNA", [ # 3 GO-headers
        "GO:0071824",   # BP **   11 uGOs    83   4 L04 D04 R04 AH    .... .... protein-DNA complex subunit organization
        "GO:0032993",   # CC **    1 uGOs    52  17 L02 D02 R02 C     .... .... protein-DNA complex
        "GO:0044815",   # CC **    3 uGOs     8   2 L02 D02 R02 C     .... .... DNA packaging complex
    ]),
    ("gene_silencing", [ # 1 GO-headers
        "GO:0016458",   # BP **    7 uGOs   114   4 L02 D07 R07 ABC   .... prd. gene silencing
    ]),
    ("cytoskeleton", [ # 5 GO-headers
        "GO:0007010",   # BP *     4 uGOs   385  12 L04 D04 R04 AH    .... .rdu cytoskeleton organization
        "GO:0005856",   # CC *     1 uGOs   273   9 L04 D05 R07 ABD   .... p... cytoskeleton
        "GO:0099080",   # CC **    3 uGOs    99   5 L01 D01 R01 K     .... .... supramolecular complex
        "GO:0099513",   # CC **    3 uGOs    43   4 L04 D05 R09 ABDEK .... .... polymeric cytoskeletal fiber
        "GO:0005198",   # MF **    1 uGOs    46  21 L01 D01 R01 G     .... .... structural molecule activity
    ]),
    ("telomere", [ # 2 GO-headers
        "GO:0032200",   # BP **    3 uGOs    65   2 L05 D05 R05 AH    .... p... telomere organization
        "GO:0000781",   # CC **    2 uGOs    10   2 L05 D06 R10 ABDE  .... p... chromosome, telomeric region
    ]),
    ("stimulus", [ # 1 GO-headers
        "GO:0050896",   # BP *    12 uGOs  6000  16 L01 D01 R01 F     .... .rdu response to stimulus
    ]),
    ("cell adhesion", [ # 2 GO-headers
        "GO:0022610",   # BP **    1 uGOs   350   4 L01 D01 R01 P     .... .... biological adhesion
        "GO:0007155",   # BP **   14 uGOs   311   8 L02 D02 R02 P     .... .rdu cell adhesion
    ]),
    ("Reproduction", [ # 2 GO-headers
        "GO:0022414",   # BP *     3 uGOs  1221 106 L01 D01 R02 LM    P... .rdu reproductive process
        "GO:1990111",   # CC **    1 uGOs     0   0 L04 D06 R06 ABC   .... .... spermatoproteasome complex
    ]),
    ("embryo", [ # 2 GO-headers
        "GO:0048856",   # BP *    13 uGOs  6065 195 L02 D02 R02 E     .... p... anatomical structure development
        "GO:0009790",   # BP *     6 uGOs   649   2 L03 D04 R04 DE    .... prdu embryo development
    ]),
    ("hormone", [ # 4 GO-headers
        "GO:0010817",   # BP **    1 uGOs   451   8 L03 D03 R03 B     .... .... regulation of hormone levels
        "GO:0009725",   # BP *     1 uGOs   271  29 L03 D04 R04 F     .... .... response to hormone
        "GO:0042445",   # BP **   12 uGOs   236  21 L02 D04 R04 BC    .... .rdu hormone metabolic process
        "GO:0005179",   # MF **    1 uGOs    19  14 L04 D05 R05 BE    .... .... hormone activity
    ]),
    ("protein_breakdown", [ # 11 GO-headers
        "GO:0006950",   # BP *     1 uGOs  1600  30 L02 D02 R02 F     .... .r.. response to stress
        "GO:0006508",   # BP **    7 uGOs   268   6 L04 D05 R05 C     .... .rdu proteolysis
        "GO:0070647",   # BP **    5 uGOs   122   2 L06 D07 R07 AC    .... .rdu protein modification by small protein conjugation or removal
        "GO:0035966",   # BP *     2 uGOs    30   3 L03 D04 R04 F     .... .... response to topologically incorrect protein
        "GO:0000502",   # CC *     2 uGOs    26   3 L03 D05 R05 ABC   .... p... proteasome complex
        "GO:0008180",   # CC **    1 uGOs     0   0 L02 D05 R09 ABCDE .... .... COP9 signalosome
        "GO:0008233",   # MF **    6 uGOs    75   4 L03 D03 R03 A     .... .... peptidase activity
        "GO:0061134",   # MF **    5 uGOs    15   3 L03 D03 R03 E     .... .... peptidase regulator activity
        "GO:0101005",   # MF **    2 uGOs     8   6 L03 D03 R03 A     .... .... ubiquitinyl hydrolase activity
        "GO:0002020",   # MF **    1 uGOs     2   2 L04 D04 R04 B     .... .... protease binding
        "GO:0051787",   # MF **    1 uGOs     0   0 L03 D03 R03 B     .... .... misfolded protein binding
    ]),
    ("peptide", [ # 5 GO-headers
        "GO:0015833",   # BP *     3 uGOs   552  12 L05 D06 R06 G     .... .r.. peptide transport
        "GO:0006518",   # BP **    2 uGOs   356  21 L04 D05 R05 AC    .... .... peptide metabolic process
        "GO:0018149",   # BP **    1 uGOs    61  34 L06 D07 R07 AC    .... .... peptide cross-linking
        "GO:0042277",   # MF **    2 uGOs    45  10 L03 D03 R03 B     .... .... peptide binding
        "GO:0015440",   # MF **    1 uGOs     8   4 L05 D12 R12 AC    .... .... peptide-transporting ATPase activity
    ]),
    ("protein", [ # 7 GO-headers
        "GO:0006464",   # BP **    8 uGOs  1319  62 L05 D06 R06 AC    .... .rdu cellular protein modification process
        "GO:0032268",   # BP **    2 uGOs   677   9 L05 D06 R06 ABC   .R.. .... regulation of cellular protein metabolic process
        "GO:0043933",   # BP **    1 uGOs   625  13 L03 D03 R03 AH    .... .... protein-containing complex subunit organization
        "GO:0065003",   # BP **    8 uGOs   455  52 L04 D04 R04 AH    .... .rdu protein-containing complex assembly
        "GO:0072376",   # BP **    2 uGOs    59   6 L02 D05 R05 CF    .... .rdu protein activation cascade
        "GO:0006457",   # BP *     3 uGOs    22   6 L02 D02 R02 A     .... .rdu protein folding
        "GO:0071823",   # BP **    3 uGOs     3   2 L04 D04 R04 AH    .... .... protein-carbohydrate complex subunit organization
    ]),
    ("carbohydrate", [ # 2 GO-headers
        "GO:0005975",   # BP **    2 uGOs   803   9 L03 D03 R03 C     .... .rdu carbohydrate metabolic process
        "GO:0071823",   # BP **    3 uGOs     3   2 L04 D04 R04 AH    .... .... protein-carbohydrate complex subunit organization
    ]),
    ("neuro", [ # 8 GO-headers
        "GO:0007399",   # BP **    3 uGOs  1271   0 L04 D04 R05 DE    .... prdu nervous system development
        "GO:0050877",   # BP *     3 uGOs   219  14 L03 D03 R03 D     .... .rdu nervous system process
        "GO:0050808",   # BP **    2 uGOs   134   7 L03 D03 R03 AH    .... prd. synapse organization
        "GO:0099601",   # BP **    1 uGOs    10   4 L04 D07 R07 ABFJ  .... .... regulation of neurotransmitter receptor activity
        "GO:1904407",   # BP **    2 uGOs     1   1 L04 D07 R07 ABC   ...U .... positive regulation of nitric oxide metabolic process
        "GO:0099602",   # MF **    2 uGOs     3   1 L03 D03 R03 E     .... .... neurotransmitter receptor regulator activity
        "GO:0033130",   # MF **    1 uGOs     0   0 L04 D04 R04 B     .... .... acetylcholine receptor binding
        "GO:0022850",   # MF **    1 uGOs     0   0 L06 D10 R10 CD    .... .... serotonin-gated cation-selective channel activity
    ]),
    ("chemokine", [ # 5 GO-headers
        "GO:0032602",   # BP **    3 uGOs    84  22 L03 D03 R03 D     .... prdu chemokine production
        "GO:0042379",   # MF **    4 uGOs    23   5 L05 D05 R05 B     .... .... chemokine receptor binding
        "GO:0019956",   # MF *     2 uGOs    10   3 L04 D04 R04 B     .... .... chemokine binding
        "GO:0001637",   # MF *     2 uGOs    10   1 L05 D06 R06 D     .... .... G-protein coupled chemoattractant receptor activity
        "GO:0042056",   # MF **    1 uGOs     1   1 L04 D05 R05 BE    .... .... chemoattractant activity
    ]),
    ("DNA_binding", [ # 1 GO-headers
        "GO:0003677",   # MF **    4 uGOs   147  21 L04 D04 R04 B     .... .... DNA binding
    ]),
    ("DNA_repair", [ # 1 GO-headers
        "GO:0006281",   # BP **    4 uGOs   135  15 L05 D07 R07 ACF   .... prdu DNA repair
    ]),
    ("nuclease", [ # 1 GO-headers
        "GO:0004518",   # MF **    4 uGOs   107   5 L04 D04 R04 A     .... .... nuclease activity
    ]),
    ("hydrolase", [ # 2 GO-headers
        "GO:0051336",   # BP **    6 uGOs   259  25 L04 D04 R04 B     .... .... regulation of hydrolase activity
        "GO:0016798",   # MF *     2 uGOs   291   3 L03 D03 R03 A     .... .... hydrolase activity, acting on glycosyl bonds
    ]),
    ("RNA", [ # 1 GO-headers
        "GO:0016070",   # BP **   16 uGOs  1224  17 L05 D06 R06 AC    .... .rdu RNA metabolic process
    ]),
    ("tRNA", [ # 3 GO-headers
        "GO:0006399",   # BP *     1 uGOs   143   7 L07 D08 R08 AC    .... .rdu tRNA metabolic process
        "GO:0006399",   # BP *     1 uGOs   143   7 L07 D08 R08 AC    .... .rdu tRNA metabolic process
        "GO:0140101",   # MF *     1 uGOs   111  41 L03 D03 R03 A     .... .... catalytic activity, acting on a tRNA
    ]),
    ("m_nit", [ # 2 GO-headers
        "GO:0006807",   # BP *     6 uGOs  6236  16 L02 D02 R02 C     .... .rdu nitrogen compound metabolic process
        "GO:0034641",   # BP *     5 uGOs  3595  23 L03 D03 R03 AC    .... .... cellular nitrogen compound metabolic process
    ]),
    ("nitric_oxide", [ # 1 GO-headers
        "GO:0018916",   # BP **    1 uGOs     2   2 L03 D05 R05 AC    .... .... nitrobenzene metabolic process
    ]),
    ("vacuole", [ # 1 GO-headers
        "GO:0005773",   # CC *     2 uGOs    85   6 L04 D05 R07 ABD   .... p... vacuole
    ]),
    ("chromosome", [ # 2 GO-headers
        "GO:0051276",   # BP **    1 uGOs   683  18 L04 D04 R04 AH    .... prdu chromosome organization
        "GO:0005694",   # CC **    5 uGOs   215   6 L04 D05 R07 ABD   .... p... chromosome
    ]),
    ("grow", [ # 3 GO-headers
        "GO:0040007",   # BP *     3 uGOs   432  11 L01 D01 R01 O     .... .rdu growth
        "GO:0010573",   # BP *     1 uGOs     3   0 L03 D03 R03 D     .... .rdu vascular endothelial growth factor production
        "GO:0070851",   # MF **    2 uGOs    30  21 L04 D04 R04 B     .... .... growth factor receptor binding
    ]),
    ("temperature", [ # 2 GO-headers
        "GO:0009266",   # BP *     2 uGOs    39   7 L03 D03 R03 F     .... .... response to temperature stimulus
        "GO:0001659",   # BP **    3 uGOs    23   5 L03 D05 R05 BD    .... .... temperature homeostasis
    ]),
    ("lysis", [ # 2 GO-headers
        "GO:0000323",   # CC *     1 uGOs    47   5 L05 D06 R08 ABD   .... p... lytic vacuole
        "GO:0005764",   # CC *     1 uGOs    33   6 L06 D07 R09 ABD   .... p... lysosome
    ]),
    ("body_fluid", [ # 1 GO-headers
        "GO:0050878",   # BP **    2 uGOs    97  10 L03 D03 R03 B     .... p... regulation of body fluid levels
    ]),
    ("pigment", [ # 2 GO-headers
        "GO:0042440",   # BP **    2 uGOs   136   9 L02 D02 R02 C     .... .... pigment metabolic process
        "GO:0033013",   # BP *     5 uGOs    95   7 L04 D04 R04 AC    .... .rdu tetrapyrrole metabolic process
    ]),
    ("glycosyl_transfer", [ # 2 GO-headers
        "GO:0016757",   # MF **    5 uGOs   557   7 L03 D03 R03 A     .... .... transferase activity, transferring glycosyl groups
        "GO:0015020",   # MF **    1 uGOs    12  12 L05 D05 R05 A     .... .... glucuronosyltransferase activity
    ]),
    ("transferase", [ # 4 GO-headers
        "GO:0016757",   # MF **    5 uGOs   557   7 L03 D03 R03 A     .... .... transferase activity, transferring glycosyl groups
        "GO:0016765",   # MF **    1 uGOs   112  83 L03 D03 R03 A     .... .... transferase activity, transferring alkyl or aryl (other than methyl) groups
        "GO:0015020",   # MF **    1 uGOs    12  12 L05 D05 R05 A     .... .... glucuronosyltransferase activity
        "GO:0004364",   # MF **    1 uGOs     0   0 L04 D04 R04 A     .... .... glutathione transferase activity
    ]),
    ("aminoglycan", [ # 1 GO-headers
        "GO:0006022",   # BP *     3 uGOs    68   5 L04 D04 R04 C     .... .... aminoglycan metabolic process
    ]),
    ("sulfer", [ # 1 GO-headers
        "GO:0006790",   # BP **    2 uGOs   419  63 L03 D03 R03 AC    .... .rdu sulfur compound metabolic process
    ]),
    ("binding", [ # 1 GO-headers
        "GO:0005488",   # MF *    18 uGOs  1880  50 L01 D01 R01 B     .... .... binding
    ]),
    ("membrane", [ # 4 GO-headers
        "GO:0043227",   # CC *     1 uGOs  1350   4 L02 D02 R02 D     .... p... membrane-bounded organelle
        "GO:0016020",   # CC **    6 uGOs  1043  15 L01 D01 R01 F     .... p... membrane
        "GO:0044425",   # CC **   10 uGOs   904  15 L01 D01 R02 FG    P... .... membrane part
        "GO:0009986",   # CC **    2 uGOs    19   0 L02 D02 R03 AB    .... p... cell surface
    ]),
    ("ATPase", [ # 3 GO-headers
        "GO:0016887",   # MF **    3 uGOs   184   2 L07 D07 R07 A     .... .... ATPase activity
        "GO:0015440",   # MF **    1 uGOs     8   4 L05 D12 R12 AC    .... .... peptide-transporting ATPase activity
        "GO:0015433",   # MF **    1 uGOs     0   0 L06 D13 R13 AC    .... .... peptide antigen-transporting ATPase activity
    ]),
    ("carboxylic", [ # 2 GO-headers
        "GO:0019752",   # BP **    2 uGOs  1307  32 L05 D05 R05 AC    .... .... carboxylic acid metabolic process
        "GO:0031406",   # MF **    6 uGOs    61  17 L04 D04 R04 B     .... .... carboxylic acid binding
    ]),
    ("phosphatase", [ # 1 GO-headers
        "GO:0016791",   # MF *     1 uGOs   179  70 L05 D05 R05 A     .... .... phosphatase activity
    ]),
    ("enzyme", [ # 1 GO-headers
        "GO:0030234",   # MF **    2 uGOs   153  25 L02 D02 R02 E     .... .... enzyme regulator activity
    ]),
    ("homeostasis", [ # 1 GO-headers
        "GO:0042592",   # BP **    6 uGOs   474   6 L03 D03 R03 B     .... p... homeostatic process
    ]),
    ("cytosol", [ # 1 GO-headers
        "GO:0005829",   # CC **    1 uGOs    86   1 L04 D04 R07 AB    .... p... cytosol
    ]),
    ("organelle", [ # 1 GO-headers
        "GO:0043226",   # CC *     2 uGOs  1928   7 L01 D01 R01 D     .... p... organelle
    ]),
    ("cofactor", [ # 1 GO-headers
        "GO:0051186",   # BP **    1 uGOs   552  24 L03 D03 R03 AC    .... .rdu cofactor metabolic process
    ]),
    ("broad", [ # 7 GO-headers
        "GO:0065007",   # BP *     4 uGOs 13067   3 L01 D01 R01 B     .... .... biological regulation
        "GO:0008152",   # BP *     7 uGOs  9804  21 L01 D01 R01 C     .... .rdu metabolic process
        "GO:0032501",   # BP *     1 uGOs  7546  68 L01 D01 R01 D     .... .rdu multicellular organismal process
        "GO:0009058",   # BP *    10 uGOs  3650  19 L02 D02 R02 C     .... .rdu biosynthetic process
        "GO:0071840",   # BP *     1 uGOs  3576   2 L01 D01 R01 H     .... .... cellular component organization or biogenesis
        "GO:0009056",   # BP *     3 uGOs  1839  10 L02 D02 R02 C     .... .rdu catabolic process
        "GO:0098772",   # MF **    1 uGOs   227   4 L01 D01 R01 E     .... .... molecular function regulator
    ]),
]
