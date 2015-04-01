

import sys

# Test local version of goatools
sys.path.insert(0, '..')

from goatools.obo_parser import GODag

def prt_paths(paths, PRT=sys.stdout):
  for path in paths:
    PRT.write('\n')
    for GO in path:
      PRT.write('{}\n'.format(GO))

def chk_results(actual_paths, expected_paths):
  for actual_path in actual_paths:
    # GOTerm -> list of Strings
    actual = [GO.id for GO in actual_path] 
    if actual not in expected_paths:
      raise Exception('ACTUAL {} NOT FOUND IN EXPECTED RESULTS\n'.format(actual))

def test_paths_to_top():
  #dag = GODag("./tests/data/mini_obo.obo")  
  dag = GODag("./data/mini_obo.obo")  
  expected_paths = [
    ['GO:0000001', 'GO:0000002', 'GO:0000005', 'GO:0000010'],
    ['GO:0000001', 'GO:0000003', 'GO:0000005', 'GO:0000010'],
    ['GO:0000001', 'GO:0000003', 'GO:0000006', 'GO:0000008', 'GO:0000010'] ]
  actual_paths = dag.paths_to_top("GO:0000010")
  chk_results(actual_paths, expected_paths)
  prt_paths(actual_paths)

if __name__ == '__main__':
  test_paths_to_top()
