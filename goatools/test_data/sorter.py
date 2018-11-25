"""Data for sorter tests."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


SECTIONS = [
    ("neuro", [
        "GO:0099531", # BP     32 L01 D01 U     presynaptic proc. in chemical synaptic Xmission
        "GO:0042391", # BP    117 L03 D03 A     regulation of membrane potential
        "GO:0050877", # BP     96 L03 D03 K     neurological system process
        "GO:0050808", # BP     20 L03 D03 CDI   synapse organization
        "GO:0007272", # BP     13 L03 D03 CD    ensheathment of neurons
        "GO:0051960", # BP    236 L04 D04 AB    regulation of nervous system development
        "GO:0050804", # BP    120 L03 D04 AB    modulation of synaptic transmission
        "GO:0097485", # BP     34 L04 D04 CD    neuron projection guidance
        "GO:0031644", # BP     30 L04 D04 AB    regulation of neurological system process
        "GO:0031175", # BP     14 L04 D04 CDI   neuron projection development
        "GO:0007399", # BP      0 L04 D04 F     nervous system development
        "GO:0050767", # BP    192 L05 D05 AB    regulation of neurogenesis
        "GO:0030182", # BP     71 L05 D05 CDF   neuron differentiation
        "GO:0099536", # BP     40 L04 D05 CDR   synaptic signaling
        "GO:0048666", # BP     29 L04 D05 CDF   neuron development
        "GO:0010001", # BP     17 L05 D05 CDF   glial cell differentiation
        "GO:0051969", # BP      5 L03 D05 AB    regulation of transmission of nerve impulse
        "GO:0022008", # BP      3 L05 D05 CDF   neurogenesis
        "GO:0014002", # BP      1 L05 D06 CDF   astrocyte development
        "GO:0048812", # BP     27 L05 D07 CDFI  neuron projection morphogenesis
        "GO:0048667", # BP      6 L06 D07 CDFI  cell morphogenesis in neuron differentiation
        "GO:0072578", # BP      5 L05 D07 CDHI  neurotransmitter-gated ion channel clustering
        "GO:0007409", # BP     23 L06 D08 CDFI  axonogenesis
    ]),
    ("immune", [
        "GO:0002376", # BP    564 L01 D01 M     immune system process
        "GO:0002682", # BP  1,183 L02 D02 AB    regulation of immune system process
        "GO:0001817", # BP    476 L03 D03 AB    regulation of cytokine production
        "GO:0001816", # BP    110 L03 D03 DK    cytokine production
        "GO:0034097", # BP     59 L04 D04 G     response to cytokine
        "GO:0045087", # BP     25 L03 D04 GM    innate immune response
        "GO:0006954", # BP     25 L04 D04 G     inflammatory response
        "GO:0002521", # BP     72 L05 D05 CDF   leukocyte differentiation
        "GO:0002467", # BP      0 L03 D05 FGM   germinal center formation
        "GO:0007229", # BP      0 L05 D05 AB    integrin-mediated signaling pathway
        "GO:0050900", # BP     57 L02 D06 CDMN  leukocyte migration
        "GO:0042130", # BP      9 L07 D08 AB    negative regulation of T cell proliferation
    ]),
    ("adhesion", [
        "GO:0022610", # BP    194 L01 D01 P     biological adhesion
        "GO:0030155", # BP    246 L02 D02 AB    regulation of cell adhesion
        "GO:0007155", # BP    165 L02 D02 P     cell adhesion
    ]),
    ("viral/bacteria", [
        "GO:0016032", # BP    301 L03 D04 CJ    viral process
        "GO:0050792", # BP    119 L03 D04 AB    regulation of viral process
        "GO:0098542", # BP     37 L03 D05 GJ    defense response to other organism
        "GO:0009617", # BP     12 L03 D05 GJ    response to bacterium
    ]),
    ("cell cycle", [
        "GO:0022402", # BP    463 L02 D02 C     cell cycle process
        "GO:0022403", # BP     46 L02 D02 S     cell cycle phase
        "GO:0051726", # BP    410 L03 D03 AB    regulation of cell cycle
        "GO:0051301", # BP     54 L03 D03 CD    cell division
        "GO:0007049", # BP     12 L03 D03 CD    cell cycle
        "GO:0007067", # BP      1 L04 D06 CI    mitotic nuclear division
    ]),
    ("cell death", [
        "GO:0010941", # BP    316 L03 D03 AB    regulation of cell death
        "GO:0008219", # BP    104 L03 D03 CD    cell death
        "GO:0060548", # BP    103 L04 D04 AB    negative regulation of cell death
        "GO:0097190", # BP     22 L04 D04 AB    apoptotic signaling pathway
        "GO:0097527", # BP      0 L04 D04 AB    necroptotic signaling pathway
        "GO:0008637", # BP      7 L05 D05 CI    apoptotic mitochondrial changes
    ]),
    ("lipid", [
        "GO:0006629", # BP    620 L03 D03 DE    lipid metabolic process
        "GO:0019216", # BP    243 L04 D04 AB    regulation of lipid metabolic process
        "GO:0032368", # BP    130 L04 D04 AB    regulation of lipid transport
        "GO:0033993", # BP    108 L04 D04 G     response to lipid
        "GO:0006869", # BP     93 L04 D05 DH    lipid transport
        "GO:0055088", # BP     10 L05 D05 A     lipid homeostasis
    ]),
]

# Consistent Increase
# a.p.p = antigen processing and presentation
USER_GOS = [
    "GO:0002376", # 564 L1 D01 M    immune system process
    "GO:0006955", # 100 L2 D02 GM   immune response
    "GO:0044406", #  22 L2 D02 JP   adhesion of symbiont to host
    "GO:0001817", # 476 L3 D03 AB   regulation of cytokine production
    "GO:0006952", #  92 L3 D03 G    defense response
    "GO:0002250", #  46 L3 D03 GM   adaptive immune response
    "GO:0019884", #   9 L3 D03 M    antigen processing and presentation of exogenous antigen
    "GO:0008284", #  87 L4 D04 AB   positive regulation of cell proliferation
    "GO:0045087", #  25 L3 D04 GM   innate immune response
    "GO:0006954", #  25 L4 D04 G    inflammatory response
    "GO:0002474", #   8 L4 D04 M    a.p.p. of peptide antigen via MHC class I
    "GO:0097527", #   0 L4 D04 AB   necroptotic signaling pathway
    "GO:0016477", # 235 L3 D05 CDN  cell migration
    "GO:0071407", # 119 L5 D05 CG   cellular response to organic cyclic compound
    "GO:0009617", #  12 L3 D05 GJ   response to bacterium
    "GO:0042590", #   2 L5 D05 M    a.p.p. of exogenous peptide antigen via MHC-I
    "GO:0032729", #   2 L5 D05 AB   positive regulation of interferon-gamma production
    "GO:0032755", #   2 L5 D05 AB   positive regulation of interleukin-6 production
    "GO:0002467", #   0 L3 D05 FGM  germinal center formation
    "GO:0032611", #   0 L5 D05 DK   interleukin-1 beta production
    "GO:0070269", #   0 L5 D05 CD   pyroptosis
    "GO:0007229", #   0 L5 D05 AB   integrin-mediated signaling pathway
    "GO:0030335", #  65 L5 D06 AB   positive regulation of cell migration
    "GO:0042742", #   8 L4 D06 GJ   defense response to bacterium
    "GO:0051607", #   3 L3 D06 CGJM defense response to virus
    "GO:0032760", #   2 L6 D06 AB   positive regulation of tumor necrosis factor production
    "GO:0035590", #   1 L6 D06 AB   purinergic nucleotide receptor signaling pathway
    "GO:0042832", #   1 L4 D06 GJ   defense response to protozoan
    "GO:0035458", #   0 L6 D06 CG   cellular response to interferon-beta
    "GO:0034123", #  18 L4 D07 AB   positive regulation of toll-like receptor signaling pathway
    "GO:0043123", #   4 L6 D07 AB   positive regulation of I-kappaB kinase/NF-kappaB signaling
    "GO:0050830", #   0 L5 D07 GJ   defense response to Gram-positive bacterium
    "GO:2001238", #   3 L6 D08 AB   pos. reg. of extrinsic apoptotic signaling pathway
    "GO:0001916", #   3 L5 D08 AB   pos. reg. of T cell mediated cytotoxicity
    "GO:0002726", #   2 L6 D08 AB   pos. reg. of T cell cytokine production
    "GO:0045651", #   1 L7 D08 AB   pos. reg. of macrophage differentiation
    "GO:0050718", #   0 L7 D10 AB   pos. reg. of interleukin-1 beta secretion
    "GO:0043280", #   7 L7 D11 AB   pos. reg. of cysteine-type endopeptidase act. in apoptosis
    "GO:0051092", #   1 L4 D11 AB   pos. reg. of NF-kappaB transcription factor activity
]


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
