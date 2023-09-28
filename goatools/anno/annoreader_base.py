"""Read an Association File and store the data in a Python object."""

import sys
import timeit
import datetime
import collections as cx
import logging

from ..anno.opts import AnnoOptions
from ..base import logger
from ..evidence_codes import EvidenceCodes
from ..godag.consts import NAMESPACE2NS
from ..gosubdag.go_tasks import get_go2parents_go2obj

__copyright__ = (
    "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
)
__author__ = "DV Klopfenstein"


# pylint: disable=useless-object-inheritance,too-many-public-methods
class AnnoReaderBase(object):
    """Reads a Gene Association File. Returns a Python object."""

    # pylint: disable=broad-except,line-too-long,too-many-instance-attributes

    tic = timeit.default_timer()

    # Expected values for a Qualifier
    exp_qualifiers = {
        # Seen in both GAF and gene2go
        "not",
        "contributes_to",
        "colocalizes_with",
    }
    valid_formats = {"gpad", "gaf", "gene2go", "id2gos"}
    exp_nss = {"BP", "MF", "CC"}

    def __init__(self, name, filename=None, **kws):
        # kws: allow_missing_symbol
        self.name = name  # name is one of valid_formats
        self.filename = filename
        self.godag = kws.get("godag")
        self.namespaces = kws.get("namespaces")
        self.evobj = EvidenceCodes()
        self.hdr = None
        self.datobj = None
        # pylint: disable=no-member
        self.associations = self._init_associations(filename, **kws)
        assert self.namespaces is None or isinstance(self.namespaces, set)

    def get_desc(self):
        """Get description"""
        return "{NAME} {NSs} {GODAG}".format(
            NAME=self.name,
            NSs="" if self.namespaces is None else ",".join(self.namespaces),
            GODAG="" if self.godag is None else "godag",
        )

    # pylint: disable=unused-argument
    def get_associations(self, taxid=None):
        """Get associations"""
        # taxid is for NCBI's gene2gos
        return self.associations

    def prt_summary_anno2ev(self, prt=sys.stdout):
        """Print annotation/evidence code summary."""
        self.evobj.prt_summary_anno2ev(self.associations, prt)

    def get_name(self):
        """Return type of annotation"""
        return self.name

    # pylint: disable=no-self-use
    def get_taxid(self):
        """Return taxid, if one was provided, otherwise return -1"""
        return -1

    # Arg, taxid, is used by NCBI's annotations, but not by gpad, gaf, etc.
    def get_ns2assc(self, taxid=None, **kws):
        """Return given associations into 3 (BP, MF, CC) dicts, id2gos"""
        return {
            ns: self._get_id2gos(nts, **kws)
            for ns, nts in self.get_ns2ntsanno().items()
        }

    # pylint: disable=unused-argument
    # Arg, taxid, is used by NCBI's annotations, but not by gpad, gaf, etc.
    def get_ns2ntsanno(self, taxid=None):
        """Split list of annotations into 3 lists: BP, MF, CC"""
        return self._get_ns2ntsanno(self.associations)

    # Used by gpad, gaf, etc., but not used by NCBI's annotation reader
    def _get_ns2ntsanno(self, annotations):
        """Split list of annotations into 3 lists: BP, MF, CC"""
        if self.name in {"gpad", "id2gos"}:
            assert (
                self.godag is not None
            ), "{T}: LOAD godag TO USE {C}::ns2ntsanno".format(
                C=self.__class__.__name__, T=self.name
            )
        ns2nts = cx.defaultdict(list)
        if self.godag:
            s_godag = self.godag
            for nta in annotations:
                if nta.GO_ID in s_godag:
                    ns2nts[nta.NS].append(nta)
        else:
            for nta in annotations:
                ns2nts[nta.NS].append(nta)
        return {ns: ns2nts[ns] for ns in ns2nts}

    def get_id2gos_nss(self, **kws):
        """Return all associations in a dict, id2gos, regardless of namespace"""
        return self._get_id2gos(self.associations, **kws)

    def get_id2gos(self, namespace=None, prt=sys.stdout, **kws):
        """Return associations from specified namespace in a dict, id2gos"""
        # pylint: disable=superfluous-parens
        if self.has_ns():  # Anno namedtuple has NS field
            nspc, assoc = self._get_1ns_assn(namespace)
            id2gos = self._get_id2gos(assoc, **kws)
            if prt:
                prt.write(f"{len(id2gos)} IDs in loaded association branch, {nspc}\n")
            return id2gos
        if prt and self.godag is None:
            logging.warning(
                "%s.get_id2gos: GODAG is None. IGNORING namespace(%s). If you are running `map_to_slim.py`, this warning can be ignored.",
                type(self).__name__,
                namespace,
            )
        id2gos = self._get_id2gos(self.associations, **kws)
        if prt:
            prt.write(f"{len(id2gos)} IDs in all associations\n")
        return id2gos

    def _get_1ns_assn(self, namespace_usr):
        """Get one namespace, given a user-provided namespace or a default"""
        # If all namespaces were loaded
        if self.namespaces is None:
            # Return user-specified namespace, if provided. Otherwise BP
            nspc = namespace_usr or self._get_biggest_namespace()
            # Return one namespace
            if nspc in set(NAMESPACE2NS.values()):
                return nspc, [nt for nt in self.associations if nt.NS == nspc]
            # Return all namespaces
            return nspc, self.associations
        # If one namespace was loaded, use that regardless of what user specfies
        if len(self.namespaces) == 1:
            nspc = next(iter(self.namespaces))
            if namespace_usr is not None and nspc != namespace_usr:
                logger.warning("IGNORING %s; ONLY %s WAS LOADED", namespace_usr, nspc)
            return nspc, self.associations
        if namespace_usr is None:
            logger.error(
                "get_id2gos: GODAG NOT LOADED. USING: %s",
                " ".join(sorted(self.namespaces)),
            )
        return namespace_usr, self.associations

    def _get_biggest_namespace(self):
        """Get the namespace with the most ontology terms"""
        nspc_ctr = cx.Counter([o.namespace for o in self.godag.values()])
        return max(nspc_ctr, key=nspc_ctr.get)

    def has_ns(self):
        """Return True if namespace field, NS exists on annotation namedtuples"""
        assert (
            self.associations
        ), f"NO ASSOCIATIONS IN file({self.filename}): {self.associations}"
        return hasattr(next(iter(self.associations)), "NS")

    def _get_id2gos(
        self,
        ntannos_usr,
        propagate_counts=False,
        relationships=None,
        prt=sys.stdout,
        **kws,
    ):
        """Return given ntannos_usr in a dict, id2gos"""
        options = AnnoOptions(self.evobj, **kws)
        # Default reduction is to remove. For all options, see goatools/anno/opts.py:
        #   * Evidence_Code == ND -> No biological data No biological Data available
        #   * Qualifiers contain NOT
        ntannos_indag = self._get_anno_in_dag(ntannos_usr)
        ntannos_m = self.reduce_annotations(ntannos_indag, options)
        dbid2goids = self.get_dbid2goids(
            ntannos_m, propagate_counts, relationships, prt
        )
        if options.b_geneid2gos:
            return dbid2goids
        return self._get_goid2dbids(dbid2goids)

    def _get_anno_in_dag(self, ntsanno):
        """Return annotations that are in the GODAG"""
        s_godag = self.godag
        return [nt for nt in ntsanno if nt.GO_ID in s_godag] if s_godag else ntsanno

    @staticmethod
    def _get_goid2dbids(dbid2goids):
        """Return dict of GO ID keys and a set of gene products as values"""
        goid2dbids = cx.defaultdict(set)
        for dbid, goids in dbid2goids.items():
            for goid in goids:
                goid2dbids[goid].add(dbid)
        return dict(goid2dbids)

    def _get_namespaces(self, nts):
        """Get the set of namespaces seen in the namedtuples."""
        return set(nt.NS for nt in nts) if self.has_ns() else set()

    # Qualifier (column 4)
    # Flags that modify the interpretation of an annotation one (or more) of NOT, contributes_to, colocalizes_with
    # This field is not mandatory;
    #     * cardinality 0, 1, >1;
    #     * for cardinality >1 use a pipe to separate entries (e.g. NOT|contributes_to)
    def prt_qualifiers(self, prt=sys.stdout):
        """Print Qualifiers: 1,462 colocalizes_with; 1,454 contributes_to; 1,157 not"""
        # 13 not colocalizes_with   (TBD: CHK - Seen in gene2go, but not gafs)
        #  4 not contributes_to     (TBD: CHK - Seen in gene2go, but not gafs)
        self._prt_qualifiers(self.associations, prt)

    @staticmethod
    def _prt_qualifiers(associations, prt=sys.stdout):
        """Print Qualifiers found in the annotations.
        QUALIFIERS:
             1,462 colocalizes_with
             1,454 contributes_to
             1,157 not
                13 not colocalizes_with   (TBD: CHK - Seen in gene2go, but not gafs)
                 4 not contributes_to     (TBD: CHK - Seen in gene2go, but not gafs)
        """
        prt.write("QUALIFIERS:\n")
        for fld, cnt in cx.Counter(
            q for nt in associations for q in nt.Qualifier
        ).most_common():
            prt.write(f"    {cnt:6,} {fld}\n")

    def reduce_annotations(self, annotations, options):
        """Reduce annotations to ones used to identify enrichment (normally exclude ND and NOT)."""
        getfnc_qual_ev = options.getfnc_qual_ev()
        return [
            nt for nt in annotations if getfnc_qual_ev(nt.Qualifier, nt.Evidence_Code)
        ]

    @staticmethod
    def update_association(assc_goidsets, go2ancestors, prt=sys.stdout):
        """Update the GO sets in assc_gene2gos to include all GO ancestors"""
        goids_avail = set(go2ancestors)
        # assc_gos is assc_gene2gos.values()
        for assc_goids_cur in assc_goidsets:
            parents = set()
            for goid in assc_goids_cur.intersection(goids_avail):
                parents.update(go2ancestors[goid])
            assc_goids_cur.update(parents)

    def _get_go2ancestors(self, goids_assoc_usr, relationships, prt=sys.stdout):
        """Return go2ancestors (set of parent GO IDs) for all GO ID keys in go2obj."""
        assert self.godag is not None
        _godag = self.godag
        # Get GO IDs in annotations that are in GO DAG
        goids_avail = set(_godag)
        goids_missing = self._rpt_goids_notfound(goids_assoc_usr, goids_avail)
        goids_assoc_cur = goids_assoc_usr & goids_avail - goids_missing
        # Get GO Term for each current GO ID in the annotations
        _go2obj_assc = {go: _godag[go] for go in goids_assoc_cur}
        go2ancestors = get_go2parents_go2obj(_go2obj_assc, relationships, prt)
        if prt:
            prt.write("{len(goids_avail)} GO IDs -> {len(go2ancestors)} go2ancestors\n")
        return go2ancestors

    @staticmethod
    def _rpt_goids_notfound(goids_assoc_all, goids_avail):
        """Report the number of GO IDs in the association, but not in the GODAG"""
        goids_missing = goids_assoc_all.difference(goids_avail)
        if goids_missing:
            print(
                "{N} GO IDs NOT FOUND IN ASSOCIATION: {GOs}".format(
                    N=len(goids_missing), GOs=" ".join(sorted(goids_missing))
                )
            )
        return goids_missing

    def get_dbid2goids(
        self, ntannos, propagate_counts=False, relationships=None, prt=sys.stdout
    ):
        """Return gene2go data for user-specified taxids."""
        if propagate_counts:
            return self._get_dbid2goids_p1(ntannos, relationships, prt)
        return self._get_dbid2goids_p0(ntannos)

    @staticmethod
    def _get_dbid2goids_p0(associations):
        """Return gene2goids with annotations as-is (propagate_counts == False)"""
        id2gos = cx.defaultdict(set)
        for ntd in associations:
            id2gos[ntd.DB_ID].add(ntd.GO_ID)
        return dict(id2gos)

    def _get_dbid2goids_p1(self, ntannos, relationships=None, prt=sys.stdout):
        """Return gene2goids with propagate_counts == True"""
        id2gos = cx.defaultdict(set)
        goids_annos = set(nt.GO_ID for nt in ntannos)
        go2ancestors = self._get_go2ancestors(goids_annos, relationships, prt)
        # https://github.com/geneontology/go-annotation/issues/3523
        ## exclude = {'GO:2000325', 'GO:2000327'}
        for ntd in ntannos:
            goid = ntd.GO_ID
            goids = id2gos[ntd.DB_ID]
            goids.add(goid)
            if goid in go2ancestors:
                goids.update(go2ancestors[goid])
        return dict(id2gos)

    @staticmethod
    def get_goid2dbids(associations):
        """Return gene2go data for user-specified taxids."""
        go2ids = cx.defaultdict(set)
        for ntd in associations:
            go2ids[ntd.GO_ID].add(ntd.DB_ID)
        return dict(go2ids)

    def hms(self, msg, tic=None, prt=sys.stdout):
        """Print elapsed time and message."""
        if tic is None:
            tic = self.tic
        now = timeit.default_timer()
        hms = str(datetime.timedelta(seconds=now - tic))
        prt.write(f"{hms}: {msg}\n")
        return now

    def chk_associations(self, fout_err=None):
        """Check that associations are in expected format."""
        # pylint: disable=unnecessary-pass
        pass

    def nts_ev_nd(self):
        """Get annotations where Evidence_code == 'ND' (No biological data)"""
        return [nt for nt in self.associations if nt.Evidence_Code == "ND"]

    def nts_qual_not(self):
        """Get annotations having Qualifiers containing NOT"""
        return [nt for nt in self.associations if self._has_not_qual(nt)]

    def chk_qualifiers(self):
        """Check format of qualifier"""
        if self.name == "id2gos":
            return
        for ntd in self.associations:
            qual = ntd.Qualifier
            assert isinstance(
                qual, set
            ), f"{self.name}: QUALIFIER MUST BE A LIST: {ntd}"
            assert qual != {""}, ntd
            assert qual != {"-"}, ntd
            assert "always" not in qual, "SPEC SAID IT WOULD BE THERE"

    def chk_godag(self):
        """Check that a GODag was loaded"""
        if not self.godag:
            raise RuntimeError(
                "{CLS} MUST INCLUDE GODag: {CLS}(file.anno, godag=godag)".format(
                    CLS=self.__class__.__name__
                )
            )

    @staticmethod
    def _has_not_qual(ntd):
        """Return True if the qualifiers contain a 'NOT'"""
        for qual in ntd.Qualifier:
            if "not" in qual:
                return True
            if "NOT" in qual:
                return True
        return False

    def prt_counts(self, prt=sys.stdout):
        """Print the number of taxids stored."""
        num_annos = len(self.associations)
        # 792,891 annotations for 3 taxids stored: 10090 7227 9606
        prt.write(f"{num_annos:8,} annotations\n")


# Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
