
import os

def test():
  """Test to find this error below.

       Traceback (most recent call last):
         File "../scripts/find_enrichment.py", line 130, in <module>
           study=study, methods=methods)
         File "../scripts/../goatools/go_enrichment.py", line 93, in __init__
           self.run_study(study)
         File "../scripts/../goatools/go_enrichment.py", line 129, in run_study
           p_val_distribution = calc_qval(study_count, study_n,
       UnboundLocalError: local variable 'study_count' referenced before assignment
  """
  os.system("python {SCR} --alpha=0.05 {STUDY} {POP} {ASSN} --fdr --obo={OBO}".format(
    SCR="../scripts/find_enrichment.py",
    OBO="../go-basic.obo",
    STUDY="data/study_unknown",
    POP="../data/population",
    ASSN="../data/association"))

if __name__ == '__main__':
  test()

