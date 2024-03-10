#!/bin/bash

#ASK FOR FAIL THRESHOLD
failtheshold=3
echo "Enter threshold for max no. ,inclusive, of fails allowed in sequencing quality report. (Input integer)"
echo "Will remove samples will more fails than threshold"
echo "(Default fail threshold = 3)"
read failthreshold

#1. Mean the replicates for each condition
##1a. create list of all unique conditions
echo "Creaing conditions list..."
touch ./temp/replicateindex.txt
awk 'BEGIN{FS="\t"; OFS="\t"}{if(NR!=1){print $2,$4,$5}}' ./seqdata/Tco2.fqfiles | 
sort - | 
uniq - > ./temp/replicateindex.txt

##1b. create reference index
echo "Creating reference index..."
touch ./temp/replicateindexref.txt
awk 'BEGIN{FS="\t";}{if(NR!=1){print $1,$2,$4,$5}}' ./seqdata/Tco2.fqfiles > ./temp/replicateindexref.txt

##1bi. Removing samples with low quality, higher fails in quality check than threshold.
echo "Removing samples with low quality sequencing in either/or end1, end2"

awk 'BEGIN{FS="\t"}{if(NR != 1){print $0}}' ./temp/reportfinal.txt > ./temp/reportfinalnocol.txt
awk -v thresh="$failthreshold" 'BEGIN{FS="\t"}{if($3 >= thresh || $6 >= thresh){print $0}}' ./temp/reportfinalnocol.txt > ./temp/failthresholdrm.txt
cut -d " " -f 1 ./temp/failthresholdrm.txt > ./temp/failthresholdrm1.txt


grep -v -f ./temp/failthresholdrm1.txt ./temp/replicateindexref.txt > ./temp/replicateindexreftemp.txt
mv -f ./temp/replicateindexreftemp.txt ./temp/replicateindexref.txt

##1c. separate search into 3 files for each condition
echo "Creating search index..."
awk 'BEGIN{FS="\t";}{print $1;}' ./temp/replicateindex.txt | sort - | uniq - > ./temp/search1.txt
awk 'BEGIN{FS="\t";}{print $2;}' ./temp/replicateindex.txt | sort - | uniq - > ./temp/search2.txt
awk 'BEGIN{FS="\t";}{print $3;}' ./temp/replicateindex.txt | sort - | uniq - > ./temp/search3.txt

##1d. group replicates together
echo "Grouping replicates..."

mkdir ./temp/search
while read search1
do
grep -w ${search1} ./temp/replicateindexref.txt > ./temp/search1temp.txt

  while read search2
  do
  grep -w ${search2} ./temp/search1temp.txt > ./temp/search2temp.txt
  
    while read search3
    do
    grep -w ${search3} ./temp/search2temp.txt > ./temp/search/${search1}_${search2}_${search3}
    echo "${search1}_${search2}_${search3} groups created"
    
    done < ./temp/search3.txt
  done < ./temp/search2.txt
done < ./temp/search1.txt

rm -f ./temp/search/*_0_Induced


##1e.gather only counts of each gene of each sample
echo "Gathering counts for each gene of each sample..."

ls ./temp/search > ./temp/searchdirindex.txt
mkdir ./counts/groups

while read file
do
echo "Gathering ${file}..."
cat ./temp/search/${file} | awk '{FS=" ";}{print $1}' - > ./temp/repgroups.txt

  while read repgroups
  do
  echo "Gathering ${repgroups}..."
  cat ./counts/${repgroups:0:3}-${repgroups:3}_counts.txt | awk '{FS="\t";}{if($6 != "hot"){print $6}}' - > ./counts/groups/${repgroups:0:3}-${repgroups:3}_counts_only.txt
  
  done < ./temp/repgroups.txt
done < ./temp/searchdirindex.txt


#1f. append counts of each condition to one file
echo "Appending counts of each condition to one file"
mkdir ./counts/reptotal

while read file
do
cat ./temp/search/${file} | awk '{FS=" ";}{print $1}' - > ./temp/repgroups.txt
>./counts/reptotal/${file}_all.txt

  while read repgroups
  do
  paste ./counts/reptotal/${file}_all.txt ./counts/groups/${repgroups:0:3}-${repgroups:3}_counts_only.txt > ./counts/reptotal/${file}_all1.txt && mv -f ./counts/reptotal/${file}_all1.txt ./counts/reptotal/${file}_all.txt

  echo "Gathering: ./counts/${repgroups:0:3}-${repgroups:3}_counts.txt"
  echo "Group: ${file}"
  echo "Into: ./counts/reptotal/${file}_all.txt"
  echo "----------------"
  
  done < ./temp/repgroups.txt
done < ./temp/searchdirindex.txt

#1g. get means of each condition
echo "Getting means of replicates of each condition..."

while read file
do
echo "Gathering ${file}..."
cat ./temp/search/${file} | awk 'BEGIN{FS=" ";}{print $1}' - > ./temp/${file}_all.txt_repgroups.txt
done < ./temp/searchdirindex.txt

ls ./counts/reptotal > ./temp/reptotalindex.txt
mkdir ./counts/replicatemeans

while read reptotalindex
do
sum=0
denom=0
>./counts/replicatemeans/${reptotalindex}

  while read gene
  do
  sum=0
  denom=0
   
    for replicate in ${gene}
    do
    sum=$(($sum + $replicate))
    denom=$((${denom} + 1))
    done
    
  mean=echo "${sum} / ${denom}" | bc -l
  echo "${sum} / ${denom}" | bc -l >> ./counts/replicatemeans/${reptotalindex}
  echo "Counts: ${gene}"
  echo "Sum: $sum"
  echo "No. of replicates: $denom"
  echo "Mean:"
  echo ${mean}
  echo "Repindex: ${reptotalindex}"
  echo "------------------------------"
  
  done < ./counts/reptotal/${reptotalindex}
done < ./temp/reptotalindex.txt


#2. Append meancounts onto .bedfile 
mkdir ./RESULTS
mkdir ./RESULTS/MeanReplicateCounts

#2a. create base .bed file
echo "Creating base .bed file..."

> ./refseqdata/base.bed 
awk 'BEGIN{FS="\t";}{print $4,$5}' ./refseqdata/TriTrypDB-46_TcongolenseIL3000_2019.bed >> ./refseqdata/base.bed

while read reptotalindex
do
echo "Creating final ${reptotalindex} mean counts file..."
paste ./refseqdata/base.bed ./counts/replicatemeans/${reptotalindex} > ./RESULTS/MeanReplicateCounts/${reptotalindex}

done < ./temp/reptotalindex.txt


#3 Onto next step 3
source ./3foldchange.sh



##quality control for fail threshold for script.sh, then remove replicate for bad reads. talk about decreased sample size in report if need to remove. 
##talk about trimming in report (basically point out --trim is there but dont use bc. papers say its bad, find papers)
#alignment scores from bowtie2 (bowtie2) >> file, then search for alignment scores and append to quality report
#fold change beginning = 0, = fold change =inf, add small amount to denom BUT talk in report how that will inflate fold change
#talk about normalization in report, even if not asked. relate to trans splicing
# run.sh file to point towards everything



