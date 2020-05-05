hapxm identifies microhaploblock variation in bam files

In working folder:
1) A bam file containing sequences mapped by 'hapx -mb'
2) A line delimited list of target sites

Requirements (in path):
1) samtools (hapxm calls samtools view/tview)
3) GNU parallel

Usage: hapxm -b bam -o out [-F exc] [-q qual]  [-u userrange] [-db -ks] -s sites
where,
bam = path to bam file of reads aligned to ref [required]
out = name of directory for output files, not a path, will be created in current directory
sites = path to file containing genomic positions to use [required]
     Provide a line delimited list of the form contigname:site-range like:
         jcf7180008454378:303-303
         jcf7180008454378::495-495
     which specifies bp 303 of the contig named "jcf7180008454378" and bps 495-495 of contig "jcf7180008454378:".
     For now, hapxm has only been tested to handle single bp "ranges" within 1 contig.
exc = integer flag value for Samtools view -F option (properties of reads to exclude) [default=2048, excludes supplementary alignments]
qual = Samtools view -q option (minimum mapping quality of included reads) [default=1, don't use 0 use >=1 instead, 0 is poorly defined]
userrange = path to a file containing user-defined microhaploblock ranges of the form
     "contig lesite resite lemh remh mhlength" where 'lesite' is left end of the target range
     in the genome, 'remh' is the right end of the microhaploblock locus, and 'mhlength' is the
     microhaploblock length (lesite and resite are unused by the algorithm so may be dummy values):
         51jcf7180007742276 11472 11472 11472 11474 3
         51jcf7180007742276 11493 11493 11491 11500 10

-db = debugging mode, save some internal data structures as files (may produce a lot of files)
-ks = keep singletons, default behavior is to ignore microhaplotype singletons (occur in only 1 sequence)

Examples: hapxm.sh -b /share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/Hs1pro1l1.finalaln.bam \
            -o hxm1 -db -s <(for i in $(seq 11423 1 11500); do echo 51jcf7180007742276:"$i"-"$i"; done;)

hapxm.sh -b /share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/Hs1pro1l1.finalaln.bam \
            -o hxm1singl -db -ks -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)

for i in $(seq 50 1 55);
  do hapxm.sh -b /share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/"$i"Hs1pro1l1.finalaln.bam.TMP \
            -o hxm"$i" -db -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)
  done;
  
for i in $(seq 50 1 55);
  do hapxm.sh -b /share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/"$i"Hs1pro1l1.finalaln.bam.TMP \
            -o hxm"$i"onhxm1 -db -u /share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/hxm1/mhendstiled.txt \
            -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)
  done;



#Postprocessing
#determine microhaploblock major alleles that differ between pools
rhead=$(grep ^"#" hxm50onhxm1/hapxmlog.txt | cut -d' ' -f1-4);
head=$(head -1 <<<"$rhead");
rhead=$(tail -n +2 <<<"$rhead");#remove col head from row heads

outmajall="";
outmajfreq="";
for i in $(seq 50 1 55);
  do majall=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f9 | cut -d: -f1); #extract the allele with highest frequency
    majfreq=$(grep ^"#" hxm"$i"onhxm1/hapxmlog.txt | tail -n +2 | cut -d' ' -f8 | cut -d: -f1); #extract the frequency of the major allele
    outmajall=$(paste -d: <(echo "$outmajall") <(echo "$majall"));
    outmajfreq=$(paste -d: <(echo "$outmajfreq") <(echo "$majfreq"));
  done;
outmajall=$(sed 's/^://' <<<"$outmajall"); #remove leading colon
outmajfreq=$(sed 's/^://' <<<"$outmajfreq"); #remove leading colon

#recode mhblocks from DNA sequence to integers
set -f; #you have to undo globbing before running this command since some values are solely asterisks
outrecode="";
while read l;
do a=$(echo "$l" | tr ':' '\n' | sort -u); #list containing DNA sequences of alleles
  #a=$(echo "$l" | tr ':' '\n' | sort -u | sed 's/^/"/' | sed 's/$/"/'); #list containing DNA sequences of alleles
  j=1; #numeric allele name
  ll="$l";
  if [[ "$a" == "" ]];
  then outrecode+="$l"$'\n'; #deal with mhloci with no alleles by placing :::::
  else
    for aa in $a; 
      do ll=${ll//"$aa"/"$j"}; #built-in bash replace obviates need to escape asterisks using sed
        #ll=$(sed 's/'"$aa"'/'$j'/g' <<<"$ll");
        j=$(( $j + 1 ));
      done;
     outrecode+="$ll"$'\n';
  fi;
done <<<"$outmajall";
outrecode=$(sed '/^$/d' <<<"$outrecode"); #remove terminal blank line
set +f; #redo globbing


#assemble output
outp1=$(paste -d' ' <(echo "$rhead") <(echo "$outmajall") <(echo "$outrecode") <(echo "$outmajfreq")); #add row header to major allele list

#find all mh loci with different major alleles, label those rows with @
outp2="";
outp2=$(while read l;
do a=$(echo "$l" | cut -d' ' -f5 | tr ':' '\n'| sed '/^$/d' | sort | uniq | wc -l);
  if [[ "$a" > 1 ]];
  then echo "@$l";
  else echo "$l"
  fi;
done <<<"$outp1";
);




