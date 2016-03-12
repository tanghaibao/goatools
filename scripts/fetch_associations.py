"""Script for downloading associations from GO Database for a given species."""
# -*- coding: UTF-8 -*-
from __future__ import print_function

"""
Fetch associations for S Pombe (NCBI Taxonomy ID 4896)
>>> python {SCR} --taxon_id 4896 -o spombe.assocs

The output format is the same as for 

    SPBC18H10.18c   GO:0008150;GO:0005737;GO:0016021;GO:0003674
    mis13   GO:0005634;GO:0005515;GO:0031617;GO:0000070;GO:0000779;GO:0000941;GO:0000444
    fap7    GO:0016887;GO:0005634;GO:0000462;GO:0004017;GO:0005524;GO:0017111;GO:0005829
    SPCC594.01      GO:0008150;GO:0005575;GO:0003674

Other useful taxon IDs:

 - 9606 Human
 - 10090 Mouse

Note that sometimes it is necessary to use the strain ID. See: http://geneontology.org/page/download-annotations

TODO: allow use of taxon labels, allow custom filtering

""".format(SCR=__file__)

import pysolr
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))



def main():
    """Fetch simple gene-term assocaitions from Golr using bioentity document type, one line per gene."""
    import argparse
    prs = argparse.ArgumentParser(__doc__,
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    prs.add_argument('--taxon_id', type=str,
                     help='NCBI taxon ID, must match exact species/strain used by GO Central, e.g. 4896 for S Pombe')
    prs.add_argument('--golr_url', default='http://golr.geneontology.org/solr/', type=str,
                     help='NCBI taxon ID, must match exact species/strain used by GO Central, e.g. 4896 for S Pombe')
    prs.add_argument('-o', default=None, type=str,
                     help="Specifies the name of the output file")
    prs.add_argument('--max_rows', default=100000, type=int,
                     help="maximum rows to be fetched")

    args = prs.parse_args()

    solr = pysolr.Solr(args.golr_url, timeout=30)

    sys.stderr.write("TAX:"+args.taxon_id+"\n")
    results = solr.search(q='document_category:"bioentity" AND taxon:"NCBITaxon:'+args.taxon_id+'"',
                          fl='bioentity_label,annotation_class_list', rows=args.max_rows)
    sys.stderr.write("NUM GENES:"+str(len(results))+"\n")
    if (len(results) ==0):
        sys.stderr.write("NO RESULTS")
        exit(1)
    if (len(results) == args.max_rows):
        sys.stderr.write("max_rows set too low")
        exit(1)

    file_out = sys.stdout if args.o is None else open(args.o, 'w')
    for r in results:
        gene_symbol = r['bioentity_label']
        sys.stderr.write(gene_symbol+"\n")
        if 'annotation_class_list' in r:
            file_out.write(r['bioentity_label']+"\t" + ';'.join(r['annotation_class_list'])+"\n")
        else:
            sys.stderr.write("no annotations for "+gene_symbol+"\n")
            

    if args.o is not None:
        file_out.close()
        sys.stdout.write("  WROTE: {}\n".format(args.o))

if __name__ == "__main__":
    main()
