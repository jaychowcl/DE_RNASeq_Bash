#!/bin/bash

#1.User defined inputs
#condition1
echo "--------------------------------------------------------------------------------------------------------------------------"

echo -e "Please use format (SampleType Time Treatment) and exact names as seen on Tco2.fqfiles\nReference Condition:"
read refsample reftime reftreat

#condition2
echo -e "Please use format (SampleType Time Treatment), only vary one field from Reference Condition, and exact names as seen on Tco2.fqfiles\nQuery Condition:"
read quersample  quertime  quertreat

#range
echo -e "Please use one exact column name from Tco2.fqfiles (SampleType OR Time OR Treatment)\nRange:"
read range

echo -e "${refsample}, ${reftime}, ${reftreat}\n${quersample}, ${quertime}, ${quertreat}\n ${range}"

echo "--------------------------------------------------------------------------------------------------------------------------"
#2. Fold change between reference and query
mkdir foldchange
>./temp/pastefoldchange.txt

#2a. create index for condi ranges
echo "Creating index for condition ranges..."

>./temp/rangecompare.txt

if test ${range} == "SampleType"
then
  cp -f ./temp/search1.txt ./temp/rangecondi.txt
  
    while read rangecondi
    do
    echo ${rangecondi}_${quertime}_${quertreat} >> ./temp/rangecompare.txt
  
    done < ./temp/rangecondi.txt
  
elif test ${range} == "Time"
then
  cp -f ./temp/search2.txt ./temp/rangecondi.txt
  
    while read rangecondi
    do
    echo ${quersample}_${rangecondi}_${quertreat} >> ./temp/rangecompare.txt
  
    done < ./temp/rangecondi.txt
  
elif test ${range} == "Treatement"
then
  cp -f ./temp/search3.txt ./temp/rangecondi.txt
  
    while read rangecondi
    do
    echo ${quersample}_${quertime}_${rangecondi} >> ./temp/rangecompare.txt
  
    done < ./temp/rangecondi.txt
else
  echo "Please input valid range."
fi



#2b. reference index made by 2a to create fold changes for range then dividing
echo "Calculating fold changes..."

while read rangecondi
do
paste ./counts/replicatemeans/${refsample}_${reftime}_${reftreat}_all.txt ./counts/replicatemeans/${rangecondi}_all.txt > ./temp/pastefoldchange.txt
>./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecondi}.txt
echo "Comparing ${rangecompare}"

  while read ref query 
  do
  echo -e "ReferenceVal:${ref} \nQueryVal:${query}"
  
  
    if test $ref -eq 0
    then
      #awk 'BEGIN{FS="\t";}{$(($2 / ($1 + 0.1)))}'>> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${quersample}_${quertime}_${quertreat}.txt
      echo "${query} / (${ref} + 1)" | bc -l >> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecondi}.txt
      echo "${rangecondi} (ref==0)"
      
    else
      #awk 'BEGIN{FS="\t";}{$(($2 / $1))}'>> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${quersample}_${quertime}_${quertreat}.txt
      echo "${query} / ${ref}" | bc -l >> ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecondi}.txt
      echo "${rangecondi} (ref!=0)"
    fi
    
  done < ./temp/pastefoldchange.txt
done < ./temp/rangecompare.txt

#3. Paste them into .bedfile to create final fold changes
echo "Pasting fold changes onto .bed file into RESULTS/FoldChanges dir"
mkdir ./RESULTS/FoldChanges

while read rangecompare
do

paste ./refseqdata/TriTrypDB-46_TcongolenseIL3000_2019.bed ./foldchange/${refsample}_${reftime}_${reftreat}_vs_${rangecompare}.txt > ./RESULTS/FoldChanges/${refsample}_${reftime}_${reftreat}_vs_${rangecompare}.txt

done < ./temp/rangecompare.txt





