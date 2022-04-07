hapxm identifies microhaploblock variation in bam files

In working folder:
1) A bam file containing sequences mapped by 'hapx -mb'
2) A line delimited list of target sites

Requirements (in path):
1) samtools (hapxm calls samtools view/tview)
2) GNU parallel

Usage: hapxm -b bam -o out [-F exc] [-q qual] [-ssh mach] [-u userrange] [-db -ks -sp -ta -va] -s sites   
where,   
bam = path to bam file of reads aligned to ref [required]   
out = name of directory for output files (not a path), will be created in current directory   
sites = path to file containing genomic positions to use [required]   
     Provide a line delimited list of the form contigname:site-range like:   
          jcf7180008454378:303-303   
          jcf7180008454378:495-495   
     which specifies bp 303 of the contig named "jcf7180008454378" and bp 495-495 of contig "jcf7180008454378:".   
     For now, hapxm has only been tested to handle single bp "ranges" within 1 contig.   
exc = integer flag value for Samtools view -F option (properties of reads to exclude) [default=2048, excludes supplementary alignments]   
qual = Samtools view -q option (minimum mapping quality of included reads) [default=1, don't use 0 use >=1 instead, 0 is poorly defined]   
mach = path to "machines" file for the gnu parallel command --sshloginfile, forces distribution across nodes   
userrange = path to a file containing user-defined microhaploblock ranges of the form   
     "contig lesite resite lemh remh mhlength" where 'lesite' is left end of the target range   
     in the reference genome, 'remh' is the right end of the microhaploblock locus, and 'mhlength' is the   
     microhaploblock length (lesite and resite are unused by the algorithm so may be dummy values,   
     see also -db):   
          51jcf7180007742276 11472 11472 11472 11474 3   
          51jcf7180007742276 11493 11493 11491 11500 10   

-db = debugging mode, save some internal data structures as files (may produce a lot of files).
     Debug output files 'mhends[var,sorted,tiled].txt' are proper format for input using -u userrange
-ks = keep singletons, default behavior is to ignore microhaplotype singletons (occur in only 1 sequence)
-sp = suppress parallel processing (sets GNU parallel --jobs=1)
-ta = calculate and then process a short tiling array across -s sites (suppressed by -u userrange)
-va = calculate and then process a variant-rich tiling array across -s sites (results written to hapxmlogvar.txt). This
     option acts after -u or -ta or default processing of all microhaplotype ranges discovered. Use with -u or -ta is
     not normal, but will work.
#  
Examples:
    hapxm.sh -b /share/space/user/patellifolia/hapxtest/hapxsummary/bwamem/Hs1pro1l1.finalaln.bam \
             -o hxm1 -db -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)
	Runs hapxm.sh in -db mode, saving output to folder named 'hxm1'. Instead of an input file for -s, uses
	process substitution <(...) syntax to generate input text like:
		 51jcf7180007742276:9889-9889
		 51jcf7180007742276:9890-9890
	
    hapxm.sh -b /share/space/user/patellifolia/hapxtest/hapxsummary/bwamem/Hs1pro1l1.finalaln.bam \
             -o hxm1singl -db -ks -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)

    for i in $(seq 50 1 55);
      do hapxm.sh -b /share/space/user/patellifolia/hapxtest/hapxsummary/bwamem/"$i"Hs1pro1l1.finalaln.bam.TMP \
                  -o hxm"$i" -db -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)
      done;

    for i in $(seq 50 1 55);
      do hapxm.sh -b /share/space/user/patellifolia/hapxtest/hapxsummary/bwamem/"$i"Hs1pro1l1.finalaln.bam.TMP \
                -o hxm"$i"onhxm1 -db -u /share/space/user/patellifolia/hapxtest/hapxsummary/bwamem/hxm1/mhendstiled.txt \
                -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)
      done;
	Runs hapxm.sh in db mode on a set of bam input files, distinguished by integer values [50-55]. Uses the "tiled"
	path through microhaplotype loci calculated by a previous run of hapxm.sh -db and saved as file mhendstiled.txt.
	The purpose is to narrow the number of microhaplotype loci considered to a manageable or biologically relevant
	fraction.
#  

#Postprocessing of output from second example above, see also https://github.com/NCGRP/mb1suppl
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
do m=$(echo "$l" | tr ':' '\n' | sed 's/^$/0/' | tr '\n' ':' | sed 's/:$//'); #substitute 0 when no allele was found
  a=$(echo "$m" | tr ':' '\n' | sort -u | grep -v 0); #list containing DNA sequences of alleles, excluding 0 (missing)
  j=1; #numeric allele name
  ll="$m";
  #if [[ "$a" == "" ]];
  #then outrecode+="$l"$'\n'; #deal with mhloci with no alleles by placing :::::
  #else
    for aa in $a; 
      do ll=${ll//"$aa"/"$j"}; #built-in bash replace obviates need to escape asterisks using sed
        j=$(( $j + 1 ));
      done;
     outrecode+="$ll"$'\n';
  #fi;
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

#write out the summary
echo "$outp2" > hxmsummary.txt;

#for fun, display differing loci horizontally (this looks like a PAUP input matrix)
echo "$outp2" | grep ^@ | cut -d' ' -f6 | \
    awk -F':' '{for (f=1;f<=NF;f++) col[f] = col[f]$f} END {for (f=1;f<=NF;f++) print col[f]}';



#create input files for the R routine MicrohaploblockHeatmap.r
c=$(echo "$outp2" | cut -d' ' -f2-8); #get mhstart mhend mhlength seqs majalleles freqs for all mh majalleles
d=$(echo "$outp2" | grep ^@ | cut -d' ' -f2-8); #get mhstart mhend mhlength seqs majalleles freqs for mh majalleles that differ between pops

for k in "$c" "$d";
  do hdr="pop mhstart mhend length majallele freq seq";
    e="";
    while read l;
      do sqc=$(echo "$l" | cut -d' ' -f4); #allele sequence column
        mac=$(echo "$l" | cut -d' ' -f5); #major alleles column
        frc=$(echo "$l" | cut -d' ' -f6); #frequencies column
        j=1; #j indexes position in horizontal allele calls and freqs
        for i in $(seq 50 1 55);
          do rhead=$(echo -n "$i ";echo "$l" | cut -d' ' -f1-3); #row header labeled by population name
            sq=$(echo "$sqc" | cut -d: -f$j); #allele sequence for population $i
            ma=$(echo "$mac" | cut -d: -f$j); #major allele for population $i
            fr=$(echo "$frc" | cut -d: -f$j); #major allele frequency for population $i
            e+="$rhead $ma $fr $sq"$'\n';
    
            j=$(( $j + 1 ));
          done;
      done<<<"$k";
    e=$(sed '/^$/d' <<<"$e");
    e="$hdr"$'\n'"$e";
    
    
    if [[ "$k" == "$c" ]];
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinAll.txt;
    elif [[ "$k" == "$d" ]]; 
    then echo "$e" | tr ' ' '\t' > hxmsummaryRinDiff.txt;
    fi;
  done;
  

#create input files for adegenet genpop object
myga() {
       i=$1; #$i is a mh span like "10981 10987"
       b=$(grep "$i" "$pd"/hxm5[0-5]onhxm1/hapxmlog.txt); #lines of mh locus from all pools
       c=$(echo "$b" | cut -d' ' -f1,7,9 | awk -F' ' '{gsub("*","Z",$3); print $0}'); #isolate read count and sequence, substitute 'Z' for '*'
       u=$(echo "$c" | cut -d' ' -f3 | tr ':' ' ' | tr ' ' '\n' | sort | uniq | tr '\n' ' '); #unique alleles for mh locus
       nu=$(echo "$u" | awk -F' ' '{print NF}'); #number of unique alleles
       an=""; #allele names
       for z in $(seq 1 1 "$nu");
         do an="$an "$(echo "$i" | tr ' ' '_')".$z";
         done;
       an=$(echo "$an" | sed 's/^ //'); #remove leading space
       
       toparallel=$(
         echo "$an";
         while read ll;
         do e=$(echo "$ll" | cut -d' ' -f2 | tr ':' '\n');
           f=$(echo "$ll" | cut -d' ' -f3 | tr ':' '\n'); #replace * with Z, asterisks are a pain
           g=$(paste -d' ' <(echo "$f") <(echo "$e")); #paste like "seq count" for each allele in this pool

           k=""
           for j in $u;
             do h=$(echo "$g" | grep ^"$j " | cut -d' ' -f2); #get count of current allele
               if [[ $h == "" ]]; then h=0; fi; #set count to zero if allele not present in pool
               k="$k $h"; #add count
             done;
           echo "$k" | sed 's/^ //';
         done <<<"$c";
       );

       echo "$toparallel" \
           | awk -F' ' '{for (f=1;f<=NF;f++) col[f] = col[f]" "$f} END {for (f=1;f<=NF;f++) print col[f]}' \
           | sed 's/^ //'; #awk transpose and report, eliminate extra space in front
 
}
export -f myga;

#gather list of microhaplotype loci
cd /share/space/user/patellifolia/hapxtest/hapxsummary/bwamem;
pd=$(pwd); export pd;
a=$(grep ^# hxm5*onhxm1/hapxmlog.txt | grep -v "mhstart mhend" | cut -d' ' -f2,3 | sort -t' ' -k1,1n | uniq);
b=$(echo "$a" | parallel --keep-order --env pd --env myga myga;)

#rotate matrix and add pool names
rh="pop"$'\n'$(seq 50 1 55); $row header
tp=$(echo "$b" | awk -F' ' '{for (f=1;f<=NF;f++) col[f] = col[f]" "$f} END {for (f=1;f<=NF;f++) print col[f]}' | sed 's/^ //');
tpm=$(paste -d' ' <(echo "$rh") <(echo "$tp"));

#write result
echo "$tpm" > AdegenetGenpopImport.txt;
