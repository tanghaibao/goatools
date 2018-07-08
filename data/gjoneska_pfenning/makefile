
CHART := david_chart6p7_Consistent_Increase.txt

cat_cnt:
	#perl -ne '@f = split(/\t/, $_); print "$f[0]\n"' $(CHART) | sort | uniq -c
	#perl -ne '@f = split(); print "$f[0]\n"' $(CHART) | sort | uniq -c
	#perl -ne '@f = split(/\t/); printf "%s\n", $f[0]' $(CHART) | sort | uniq -c
	perl -ne '@f = split(/\t/); printf "%s\n", @f' $(CHART) | sort | uniq -c
