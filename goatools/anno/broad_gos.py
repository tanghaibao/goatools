"""List of the GO IDs that have lots of descendants and low information content"""

__copyright__ = "Copyright (C) 2018-2019, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=line-too-long
NS2GOS_SHORT = {
    'BP': {
        'GO:0008150',  # BP 29685  18,453 0.015902  4.14 L00 D00       biological_process
        'GO:0065007',  # BP 12809  12,729 0.010970  4.51 L01 D01 A     biological regulation
        'GO:0050789',  # BP 11559  12,053 0.010387  4.57 L02 D02 A     regulation of biological process
        'GO:0009987',  # BP 11246  15,270 0.013159  4.33 L01 D01 B     cellular process
        'GO:0050794',  # BP  8212  11,011 0.009489  4.66 L03 D03 A     regulation of cellular process
        'GO:0008152',  # BP  6394   8,835 0.007614  4.88 L01 D01 C     metabolic process
        'GO:0071704',  # BP  6051   8,306 0.007158  4.94 L02 D02 C     organic substance metabolic process
        'GO:0044237',  # BP  5691   7,989 0.006885  4.98 L02 D02 BC    cellular metabolic process
    },
    'CC':{
        'GO:0005575',  # CC  4197  19,578 0.059148  2.83 L00 D00       cellular_component
        'GO:0044464',  # CC  3298  17,521 0.052934  2.94 L01 D01 A     cell part
        'GO:0044424',  # CC  2356  15,063 0.045508  3.09 L02 D02 A     intracellular part
    },
    'MF':{
        'GO:0003674',  # MF 11120  17,426 0.068263  2.68 L00 D00       molecular_function
        'GO:0003824',  # MF  7659   5,710 0.022368  3.80 L01 D01 A     catalytic activity
        'GO:0005488',  # MF  1887  15,267 0.059805  2.82 L01 D01 B     binding
    },
}

NS2GOS = {
    'BP': {
        'GO:0008150',  # BP 29685  18,453 0.015902  4.14 L00 D00       biological_process
        'GO:0065007',  # BP 12809  12,729 0.010970  4.51 L01 D01 A     biological regulation
        'GO:0050789',  # BP 11559  12,053 0.010387  4.57 L02 D02 A     regulation of biological process
        'GO:0009987',  # BP 11246  15,270 0.013159  4.33 L01 D01 B     cellular process
        'GO:0050794',  # BP  8212  11,011 0.009489  4.66 L03 D03 A     regulation of cellular process
        'GO:0008152',  # BP  6394   8,835 0.007614  4.88 L01 D01 C     metabolic process
        'GO:0071704',  # BP  6051   8,306 0.007158  4.94 L02 D02 C     organic substance metabolic process
        'GO:0044237',  # BP  5691   7,989 0.006885  4.98 L02 D02 BC    cellular metabolic process
        'GO:0044238',  # BP  4213   7,935 0.006838  4.99 L02 D02 C     primary metabolic process
        'GO:0006807',  # BP  3954   7,465 0.006433  5.05 L02 D02 C     nitrogen compound metabolic process
        'GO:0048518',  # BP  3575   6,212 0.005353  5.23 L03 D03 A     positive regulation of biological process
        'GO:0048519',  # BP  3483   5,751 0.004956  5.31 L03 D03 A     negative regulation of biological process
        'GO:0019222',  # BP  3356   7,289 0.006282  5.07 L03 D03 A     regulation of metabolic process
        'GO:0032502',  # BP  3217   4,946 0.004262  5.46 L01 D01 D     developmental process
        'GO:0031323',  # BP  2921   6,302 0.005431  5.22 L04 D04 A     regulation of cellular metabolic process
        'GO:1901564',  # BP  2881   5,283 0.004553  5.39 L03 D03 C     organonitrogen compound metabolic process
        'GO:0051239',  # BP  2674   3,248 0.002799  5.88 L03 D03 A     regulation of multicellular organismal process
        'GO:1901360',  # BP  2602   3,753 0.003234  5.73 L03 D03 C     organic cyclic compound metabolic process
        'GO:0080090',  # BP  2477   6,230 0.005369  5.23 L04 D04 A     regulation of primary metabolic process
        'GO:0048523',  # BP  2396   4,843 0.004174  5.48 L04 D04 A     negative regulation of cellular process
        'GO:0043170',  # BP  2385   6,619 0.005704  5.17 L03 D03 C     macromolecule metabolic process
        'GO:0048522',  # BP  2378   5,479 0.004722  5.36 L04 D04 A     positive regulation of cellular process
        'GO:0048583',  # BP  2372   4,334 0.003735  5.59 L03 D03 A     regulation of response to stimulus
        'GO:0044281',  # BP  2364   1,744 0.001503  6.50 L02 D02 C     small molecule metabolic process
        'GO:0050896',  # BP  2283   5,768 0.004971  5.30 L01 D01 E     response to stimulus
        'GO:0051171',  # BP  2234   6,060 0.005222  5.25 L04 D04 A     regulation of nitrogen compound metabolic process
        'GO:0006725',  # BP  2233   3,528 0.003040  5.80 L03 D03 BC    cellular aromatic compound metabolic process
        'GO:0034641',  # BP  2155   3,827 0.003298  5.71 L03 D03 BC    cellular nitrogen compound metabolic process
        'GO:0046483',  # BP  2152   3,483 0.003002  5.81 L03 D03 BC    heterocycle metabolic process
        'GO:0051179',  # BP  2119   4,811 0.004146  5.49 L01 D01 F     localization
        'GO:0060255',  # BP  2110   6,771 0.005835  5.14 L04 D04 A     regulation of macromolecule metabolic process
        'GO:0050793',  # BP  2040   2,693 0.002321  6.07 L03 D03 A     regulation of developmental process
        'GO:0009058',  # BP  1843   2,533 0.002183  6.13 L02 D02 C     biosynthetic process
        'GO:0032879',  # BP  1798   2,851 0.002457  6.01 L03 D03 A     regulation of localization
        'GO:1901576',  # BP  1766   2,475 0.002133  6.15 L03 D03 C     organic substance biosynthetic process
        'GO:0071840',  # BP  1738   5,385 0.004641  5.37 L01 D01 G     cellular component organization or biogenesis
        'GO:0044260',  # BP  1726   4,739 0.004084  5.50 L03 D04 BC    cellular macromolecule metabolic process
        'GO:0065008',  # BP  1710   3,759 0.003239  5.73 L02 D02 A     regulation of biological quality
        'GO:0016043',  # BP  1709   5,337 0.004599  5.38 L02 D02 BG    cellular component organization
        'GO:0051234',  # BP  1681   4,325 0.003727  5.59 L02 D02 F     establishment of localization
        'GO:0006810',  # BP  1587   4,195 0.003615  5.62 L03 D03 F     transport
        'GO:0009889',  # BP  1583   4,400 0.003792  5.57 L04 D04 A     regulation of biosynthetic process
        'GO:0044249',  # BP  1575   2,346 0.002022  6.20 L03 D03 BC    cellular biosynthetic process
        'GO:2000026',  # BP  1532   2,154 0.001856  6.29 L04 D04 A     regulation of multicellular organismal development
        'GO:0006139',  # BP  1524   3,312 0.002854  5.86 L03 D04 BC    nucleobase-containing compound metabolic process
        'GO:0051704',  # BP  1474   1,534 0.001322  6.63 L01 D01 H     multi-organism process
        'GO:0031326',  # BP  1467   4,321 0.003724  5.59 L05 D05 A     regulation of cellular biosynthetic process
        'GO:0009056',  # BP  1382   1,946 0.001677  6.39 L02 D02 C     catabolic process
        'GO:0051128',  # BP  1372   2,507 0.002160  6.14 L04 D04 A     regulation of cellular component organization
        'GO:0010646',  # BP  1323   3,537 0.003048  5.79 L04 D04 A     regulation of cell communication
        'GO:0023051',  # BP  1302   3,579 0.003084  5.78 L03 D03 A     regulation of signaling
        'GO:1901575',  # BP  1290   1,678 0.001446  6.54 L03 D03 C     organic substance catabolic process
        'GO:0051049',  # BP  1280   1,879 0.001619  6.43 L04 D04 A     regulation of transport
        'GO:0019538',  # BP  1171   4,189 0.003610  5.62 L03 D04 C     protein metabolic process
        'GO:0044248',  # BP  1164   1,693 0.001459  6.53 L03 D03 BC    cellular catabolic process
        'GO:0006082',  # BP  1130     998 0.000860  7.06 L03 D03 BC    organic acid metabolic process
        'GO:0002682',  # BP  1127   1,659 0.001430  6.55 L03 D03 A     regulation of immune system process
        'GO:0043436',  # BP  1089     977 0.000842  7.08 L04 D04 BC    oxoacid metabolic process
        'GO:0048869',  # BP  1062   2,764 0.002382  6.04 L02 D02 BD    cellular developmental process
        'GO:0042221',  # BP  1059   2,723 0.002347  6.05 L02 D02 E     response to chemical
        'GO:0009893',  # BP  1056   3,640 0.003137  5.76 L04 D04 A     positive regulation of metabolic process
        'GO:0043412',  # BP  1048   3,230 0.002784  5.88 L04 D04 C     macromolecule modification
        'GO:0048856',  # BP  1038   3,320 0.002861  5.86 L02 D02 D     anatomical structure development
        'GO:0009892',  # BP  1035   3,414 0.002942  5.83 L04 D04 A     negative regulation of metabolic process
        'GO:0065009',  # BP  1028   2,987 0.002574  5.96 L02 D02 A     regulation of molecular function
        'GO:0071702',  # BP  1017   2,082 0.001794  6.32 L04 D04 F     organic substance transport
        'GO:0009966',  # BP  1001   3,163 0.002726  5.90 L04 D05 A     regulation of signal transduction
        'GO:0032501',  # BP   995   3,479 0.002998  5.81 L01 D01 I     multicellular organismal process
        'GO:0019219',  # BP   987   4,138 0.003566  5.64 L05 D05 A     regulation of nucleobase-containing compound metabolic process
        'GO:0010556',  # BP   970   4,134 0.003563  5.64 L05 D05 A     regulation of macromolecule biosynthetic process
        'GO:0044267',  # BP   936   3,376 0.002909  5.84 L04 D05 BC    cellular protein metabolic process
        'GO:0031325',  # BP   936   3,328 0.002868  5.85 L05 D05 A     positive regulation of cellular metabolic process
        'GO:0051246',  # BP   936   2,935 0.002529  5.98 L05 D05 A     regulation of protein metabolic process
        'GO:0045595',  # BP   921   1,876 0.001617  6.43 L04 D04 A     regulation of cell differentiation
        'GO:0036211',  # BP   906   3,021 0.002603  5.95 L04 D05 C     protein modification process
        'GO:0006464',  # BP   904   3,021 0.002603  5.95 L05 D06 BC    cellular protein modification process
        'GO:0006793',  # BP   891   2,053 0.001769  6.34 L03 D03 BC    phosphorus metabolic process
        'GO:1901362',  # BP   875   1,318 0.001136  6.78 L04 D04 C     organic cyclic compound biosynthetic process
        'GO:0031324',  # BP   873   2,635 0.002271  6.09 L05 D05 A     negative regulation of cellular metabolic process
        'GO:0080134',  # BP   855   1,644 0.001417  6.56 L04 D04 A     regulation of response to stress
        'GO:0010468',  # BP   846   5,112 0.004405  5.42 L05 D05 A     regulation of gene expression
        'GO:0051716',  # BP   846   2,935 0.002529  5.98 L02 D02 BE    cellular response to stimulus
        'GO:0022414',  # BP   843   1,389 0.001197  6.73 L01 D01 J     reproductive process
        'GO:0090304',  # BP   833   2,754 0.002373  6.04 L04 D05 BC    nucleic acid metabolic process
        'GO:0006796',  # BP   828   2,026 0.001746  6.35 L04 D04 BC    phosphate-containing compound metabolic process
        'GO:0051240',  # BP   828   1,797 0.001549  6.47 L04 D04 A     positive regulation of multicellular organismal process
        'GO:0048584',  # BP   824   2,441 0.002104  6.16 L04 D04 A     positive regulation of response to stimulus
    },
    'CC':{
        'GO:0005575',  # CC  4197  19,578 0.059148  2.83 L00 D00       cellular_component
        'GO:0044464',  # CC  3298  17,521 0.052934  2.94 L01 D01 A     cell part
        'GO:0044424',  # CC  2356  15,063 0.045508  3.09 L02 D02 A     intracellular part
        'GO:0032991',  # CC  2107   6,150 0.018580  3.99 L01 D01 B     protein-containing complex
        'GO:0044422',  # CC  1610  10,111 0.030547  3.49 L01 D01 C     organelle part
        'GO:0044446',  # CC  1458   9,810 0.029638  3.52 L02 D03 AC    intracellular organelle part
        'GO:0044444',  # CC  1257   9,735 0.029411  3.53 L03 D03 A     cytoplasmic part
        'GO:0044425',  # CC   908   6,938 0.020961  3.87 L01 D01 D     membrane part
        'GO:0044428',  # CC   490   5,022 0.015172  4.19 L03 D04 AC    nuclear part
        'GO:0044459',  # CC   472   2,985 0.009018  4.71 L02 D02 AD    plasma membrane part
        'GO:0043226',  # CC   383  11,684 0.035299  3.34 L01 D01 E     organelle
    },
    'MF':{
        'GO:0003674',  # MF 11120  17,426 0.068263  2.68 L00 D00       molecular_function
        'GO:0003824',  # MF  7659   5,710 0.022368  3.80 L01 D01 A     catalytic activity
        'GO:0016740',  # MF  2457   2,272 0.008900  4.72 L02 D02 A     transferase activity
        'GO:0016491',  # MF  2368     748 0.002930  5.83 L02 D02 A     oxidoreductase activity
        'GO:0005488',  # MF  1887  15,267 0.059805  2.82 L01 D01 B     binding
        'GO:0016787',  # MF  1643   2,533 0.009923  4.61 L02 D02 A     hydrolase activity
        'GO:0005215',  # MF  1079   1,140 0.004466  5.41 L01 D01 C     transporter activity
        'GO:0022857',  # MF  1042   1,045 0.004094  5.50 L02 D02 C     transmembrane transporter activity
        'GO:0005515',  # MF   966  11,807 0.046252  3.07 L02 D02 B     protein binding
        'GO:0015075',  # MF   691     880 0.003447  5.67 L03 D03 C     ion transmembrane transporter activity
        'GO:0015318',  # MF   635     817 0.003200  5.74 L03 D03 C     inorganic molecular entity transmembrane transporter activity
        'GO:0016772',  # MF   603     918 0.003596  5.63 L03 D03 A     transferase activity, transferring phosphorus-containing groups
        'GO:0016788',  # MF   558     741 0.002903  5.84 L03 D03 A     hydrolase activity, acting on ester bonds
        'GO:0097159',  # MF   497   6,328 0.024789  3.70 L02 D02 B     organic cyclic compound binding
        'GO:1901363',  # MF   457   6,241 0.024448  3.71 L02 D02 B     heterocyclic compound binding
        'GO:0140096',  # MF   443   2,165 0.008481  4.77 L02 D02 A     catalytic activity, acting on a protein
        'GO:0005102',  # MF   431   1,624 0.006362  5.06 L03 D03 B     signaling receptor binding
        'GO:0060089',  # MF   419   1,538 0.006025  5.11 L01 D01 D     molecular transducer activity
        'GO:0038023',  # MF   403   1,484 0.005813  5.15 L02 D02 D     signaling receptor activity
        'GO:0016301',  # MF   355     762 0.002985  5.81 L04 D04 A     kinase activity
    },
}

# Copyright (C) 2018-2019, DV Klopfenstein. All rights reserved
