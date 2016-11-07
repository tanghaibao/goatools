import collections
import mygene


GO_KEYS_FULL = 'go.BP', 'go.MF', 'go.CC'
GO_KEYS_SPLIT = [x.split('.')[1] for x in GO_KEYS_FULL]


def parse_mygene_output(mygene_output):
    """Convert mygene.querymany output to a gene id to go term mapping (dictionary)

    Parameters
    ----------
    mygene_output : dict or list
        Dictionary (returnall=True) or list (returnall=False) of
        output from mygene.querymany

    Output
    ------
    gene_name_to_go : dict
        Mapping of gene name to a set of GO ids
    """
    # if "returnall=True" was specified, need to get just the "out" key
    if isinstance(mygene_output, dict):
        mygene_output = mygene_output['out']

    gene_name_to_go = collections.defaultdict(set)

    for line in mygene_output:
        gene_name = line['query']
        try:
            go_output = line['go']
        except KeyError:
            continue
        for go_key in GO_KEYS_SPLIT:
            try:
                go_terms = go_output[go_key]
            except KeyError:
                continue
            if isinstance(go_terms, dict):
                go_ids = set([go_terms['id']])
            else:
                go_ids = set(x['id'] for x in go_terms)
        gene_name_to_go[gene_name] |= go_ids
    return gene_name_to_go


def gene_ids_to_go(gene_ids, species='human,mouse,rat',
                   scopes='entrezgene,ensemblgene,retired,symbol',
                   fields=GO_KEYS_FULL,
                   **kwargs):
    """Get associated GO terms for each gene ID

    gene_ids : iterable of ids
        List of gene ids that you want to map
    species : str
        Comma-separated species to limit search. Default is "human,mouse,rat"
    scopes : str
        Comma-separated type of gene ids that you are giving.
        Default is "entrezgene,ensemblgene,retired,symbol"
    fields : iterable
        GO terms to use. Default is ['go.BP', 'go.MF', 'go.CC']

    Returns
    -------
    gene_to_go : dict
        Mapping of each provided gene id to a set object of GO terms
    """

    mg = mygene.MyGeneInfo()

    mygene_output = mg.querymany(gene_ids, fields=fields, scopes=scopes,
                                 species=species, **kwargs)

    gene_name_to_go = parse_mygene_output(mygene_output)
    return gene_name_to_go
