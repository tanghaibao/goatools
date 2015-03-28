
import os

def test():
  """Ensure that a study with only unknown GO Terms will run gracefully."""
  os.system("python {SCR} --alpha=0.05 {STUDY} {POP} {ASSN} --fdr --obo={OBO}".format(
    SCR="../scripts/find_enrichment.py",
    OBO="../go-basic.obo",
    STUDY="data/study_unknown",
    POP="../data/population",
    ASSN="../data/association"))

if __name__ == '__main__':
  test()

