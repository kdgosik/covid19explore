# covid19explore


## Data Splits

### NSAIDS Splits

cat pubmed_result_nsaids.txt | sed -n 1,1400044p > pubmed_result_nsaids1.txt
cat pubmed_result_nsaids.txt | sed -n 1400046,2815771p > pubmed_result_nsaids2.txt

### Hypertension Splits

cat pubmed_result_hypertension.txt | sed -n 1,1400010p > pubmed_result_hypertension1.txt
cat pubmed_result_hypertension.txt | sed -n 1400012,2800024p > pubmed_result_hypertension2.txt
cat pubmed_result_hypertension.txt | sed -n 2800026,4200025p > pubmed_result_hypertension3.txt
cat pubmed_result_hypertension.txt | sed -n 4200027,5700027p > pubmed_result_hypertension4.txt
cat pubmed_result_hypertension.txt | sed -n 5700029,7622416p > pubmed_result_hypertension5.txt