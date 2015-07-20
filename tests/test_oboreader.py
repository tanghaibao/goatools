"""Tests alternate OBOReader."""

import sys
import timeit
import datetime

# Test local version of goatools
sys.path.insert(0, '..')
from goatools.obo_parser import GODag

#################################################################
# Sub-routines to tests
#################################################################
def test_write_hier_all(dag, out):
    """Prints the entire mini GO hierarchy, with counts of children."""
    out.write('\nTEST ALL: Print all hierarchies:\n')
    dag.write_hier("GO:0000001", out, num_child=True)


def test_write_hier_norep(dag, out):
    """Shortens hierarchy report by only printing branches once.

         Prints the 'entire hierarchy' of GO:0000005 the 1st time seen:

           ---     1 GO:0000005    L-02    D-02
           ----     0 GO:0000010   L-03    D-04

         Prints just GO:0000005 (ommit child GO:10) the 2nd time seen:

           ===     1 GO:0000005    L-02    D-02

         '=' is used in hierarchy mark to indicate that the pathes
             below the marked term have already been printed.
    """
    out.write('\nTEST ALL: Print branches just once:\n')
    dag.write_hier("GO:0000001", out, num_child=True, short_prt=True)


def test_write_hier_lim(dag, out):
    """Limits hierarchy list to GO Terms specified by user."""
    go_omit = ['GO:0000005', 'GO:0000010']
    go_ids = [go_id for go_id in dag if go_id not in go_omit]
    out.write('\nTEST OMIT: 05->10:\n')
    dag.write_hier("GO:0000001", out, include_only=go_ids)
      #go_marks=[oGO.id for oGO in oGOs_in_cluster])

def test_write_hier_mrk(dag, out):
    """Print all paths, but mark GO Terms of interest. """
    mark_lst = ['GO:0000001', 'GO:0000003', 'GO:0000006', 'GO:0000008', 'GO:0000009']
    out.write('\nTEST MARK: 01->03->06->08->09:\n')
    dag.write_hier("GO:0000001", out, go_marks=mark_lst)
      #go_marks=[oGO.id for oGO in oGOs_in_cluster])

def test_all(dag_fin, use_alt, fout):
    """Run numerous tests for various reports."""
    tic = timeit.default_timer()
    dag = GODag(dag_fin, OBOReader_test=use_alt)
    toc = timeit.default_timer()
    out = sys.stdout if fout is None else open(fout, 'w')
    test_write_hier_all(dag, out)
    test_write_hier_norep(dag, out)
    test_write_hier_lim(dag, out)
    test_write_hier_mrk(dag, out)
    msg = "Elapsed HMS: {}\n\n".format(str(datetime.timedelta(seconds=(toc-tic))))
    if fout is not None:
        out.close()
        sys.stdout.write("  WROTE: {}\n".format(fout))
    sys.stdout.write(msg)

def test_run(msg, tag, use_alt=False):
    """Test hierarchy report using original/alternate OBOReader."""
    sys.stdout.write("\n{MSG}\n\n".format(MSG=msg))

    # Run test on mini-obo for quick testing
    test_all("./data/mini_obo.obo", use_alt,
             "OBOReader_{TAG}_miniobo.tmp".format(TAG=tag))


def test_oboreader_equal_(dag_fin):
    """Test that the contents of the original DAG and the alternate DAG are the same."""
    sys.stdout.write("\n\nTEST GODag EQUALITY USING {} ...\n\n".format(dag_fin))
    # 1. Read the obo file using the alternate  OBOReader.
    tic = timeit.default_timer()
    dag_alt = GODag(dag_fin, True)
    sys.stdout.write("Alternate OBOReader Elapsed HMS: {}\n\n".format(
        str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
    # 2. Read the obo file using the original OBOReader.
    tic = timeit.default_timer()
    dag_orig = GODag(dag_fin)
    sys.stdout.write("Original OBOReader Elapsed HMS: {}\n\n".format(
        str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
    # 3. Test that the contents of each GODag are the same.
    if len(dag_orig) != len(dag_alt): raise Exception("LENGTHES NOT THE SAME")
    for goid, goorig in dag_orig.items():
        goalt = dag_alt[goid]
        if goorig.id != goalt.id: raise Exception("id MISMATCH")
        if goorig.name != goalt.name: raise Exception("name MISMATCH")
        if goorig.namespace != goalt.namespace: raise Exception("namespace MISMATCH")
        if goorig.level != goalt.level: raise Exception("level MISMATCH")
        if goorig.depth != goalt.depth: raise Exception("depth MISMATCH")
        if goorig.is_obsolete != goalt.is_obsolete: raise Exception("is_obsolete MISMATCH")

        # 3a. Check that lengths of arrays are equal
        if len(goorig.parents)  != len(goalt.parents):  raise Exception("parents len MISMATCH")
        if len(goorig.children) != len(goalt.children): raise Exception("children len MISMATCH")
        if len(goorig.alt_ids)  != len(goalt.alt_ids):  raise Exception("alt_ids len MISMATCH")

        # 3b. Check that id or value of elements in array are equal
        if sorted(set(obj.id for obj in goorig.parents))  != sorted(set(obj.id for obj in goalt.parents)): 
          raise Exception("parents MISMATCH")
        if sorted(set(obj.id for obj in goorig.children)) != sorted(set(obj.id for obj in goalt.children)): 
          raise Exception("children MISMATCH")
        if sorted(set(goorig.alt_ids))  != sorted(set(goalt.alt_ids)): 
          raise Exception("alt_ids MISMATCH")

    # 4. Compare GO id keys 
    goids_orig = set(dag_orig.keys())
    goids_alt = set(dag_alt.keys())
    goids_diff = goids_alt.difference(goids_orig)
    if len(goids_diff) != 0:
        raise Exception("GODag GO id key DIFFERENCES: {}\n".format(goids_alt.difference(goids_orig)))

    # 5. Test passes
    sys.stdout.write("TEST PASSES({}): Original GODag and Alternate GODag are equal.\n".format(
      dag_fin))




#################################################################
# Tests
#################################################################
def test_oboreader_hier_rpt():
    """Compare original OBOReader with alternate OBOReader."""
    test_run("TESTING ORIGINAL OBOReader...", "orig")
    test_run("TESTING ALTERNATE OBOReader...", "alt", True)

def test_oboreader_equal():
    """Tests that contents of original DAG and alternate DAG are the same."""
    test_oboreader_equal_("./data/mini_obo.obo")
    test_oboreader_equal_("./go-basic.obo")


#################################################################
# main
#################################################################
if __name__ == '__main__':
    test_oboreader_hier_rpt()
    test_oboreader_equal()


