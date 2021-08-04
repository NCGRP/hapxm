#!/bin/bash

# Usage: see README.txt

##SUBROUTINES##
#myparsecigar() returns the distance from the contig:site-range right end to the each mapped
#paired read right end
#myparsecigar() {
#               j=$1; #pos:cigar like 11423:4M1D48M2D7M1D38M1D38M72S
#               pos=$(cut -d: -f1 <<<"$j"); #position of first aligned base
#               cig=$(cut -d: -f2 <<<"$j"); #cigar string
#               cigops=$(sed 's/[0-9]*//g' <<<"$cig" | grep -o . | sort -u | tr -d '\n'); #gather all cigar string 'operators'
#               if [[ "$cigops" == *"P"* ]];
#               then echo "CIGAR string $cig contains invalid character P, skipping";
#                 return; #skip the current contig:site-range
#               else sed 's:\([DIMSX=]\):\1\n:g' <<<"$cig" | sed '/^$/d' | grep -v "I" | sed 's/[DIMSX=]//' | awk '{s+=$1}END{print s}';
#               fi;
#}
#export -f myparsecigar;
#
#mygetends() calculates the span of bases around a contig:site-range that can be treated as a microhaplotype
#this amounts to the sequence span that is alignable across (i.e. "common to") all paired read
#haplotypes extending left and right from the contig:site-range
mygetends() {
              i=$1; #contig:site-range
              contig=$(cut -d: -f1 <<<"$i"); #contig name
              lesr=$(cut -d: -f2 <<<"$i" | cut -d'-' -f1); #left end site range
              resr=$(cut -d: -f2 <<<"$i" | cut -d'-' -f2); #right end site range
              
              #f=$(samtools view -q 1 Hs1pro1l1.finalaln.bam 51jcf7180007742276:"$i"-"$i" | cut -d$'\t' -f4 | sort -nr | head -1);
              
              #get LE and cigar string for each sequence
              g=$(samtools view -F "$stF" -q "$stq" "$bam" "$i" | cut -d$'\t' -f4,6); #-F 2048 excludes supplementary alignments
              #g=$(/share/apps/samtools view -F "$stF" -q "$stq" "$bam" "$i" | cut -d$'\t' -f4,6); #-F 2048 excludes supplementary alignments
              if [[ "$g" == "" ]]; then return; fi; #bail out if there are no aligned reads at the position
              
              #calculate closest ends to left side of site range
              le=$(cut -d$'\t' -f1 <<<"$g" | sort -nr | head -1); #position of left read pair end nearest to left end of site range
              led=$(( $lesr - $le )); #distance to closest LE
              
              #parse cigar string for right end positions
              #allres=$(echo "$g" | tr '\t' ':' | parallel --env myparsecigar myparsecigar); #variable to hold all right end positions
              allres="";
              for j in $(echo "$g" | tr '\t' ':' | tr '\n' ' ');
                do lpos=$(cut -d: -f1 <<<"$j"); #position of first aligned base
                  cig=$(cut -d: -f2 <<<"$j"); #cigar string
                  cigops=$(sed 's/[0-9]*//g' <<<"$cig" | grep -o . | sort -u | tr -d '\n'); #gather all cigar string 'operators'
                  
                  #process the cigar string. 1) split string on operations DIMSX=, 2) remove trailing blank line
                  #3) remove any line with I (insertion) op, you are looking for the right end position relative to
                  #the reference so insertions in the read are not counted, 4) remove DIMSX= characters, leaving the number
                  #of base pairs for each op, 5) count all lines, each containing the number of base pairs per op,
                  #6) subtract 1 because the cigar ops start on the lpos so that rpos is the last actual base of the 
                  #read.
                  if [[ "$cigops" == *"P"* ]];
                  then echo "CIGAR string $cig contains invalid character P, skipping";
                    return; #skip the current contig:site-range
                  else lengthcig=$(sed 's:\([DIMSX=]\):\1\n:g' <<<"$cig" | sed '/^$/d' | grep -v "I" | sed 's/[DIMSX=]//' | awk '{s+=$1}END{print s}');
                    rpos=$(( $lpos + $lengthcig - 1 )); #position of read pair right end
                    allres+="$rpos"$'\n';
                   fi;
                done;
              allres=$(sed '/^$/d' <<<"$allres"); #remove trailing blank line
              
              #calculate closest ends to right side of site range
              #distance to all REs from RE of site range
              
              re=$(sort -n <<<"$allres" | head -1); #position of closest RE
              red=$(( $re - $resr )); #distance to closest RE
              
              #report microhaplotype range to calling statement
              echo "$contig:$lesr-$resr $le-$re "$(( $re - $le + 1 ));
}
export -f mygetends;

#extract longest haploblocks across the tiling array
#for microhaploblock ranges with the same start point, keep the longest one (which is the last with current sort order)
mymhends2() {
            i=$1;
            awk -F' ' -v i=$i '$4==i{print $0}' "$pd"/mhends1.tmp | sort -t' ' -k5,5n | tail -1;
}
export -f mymhends2;

#for microhaploblock ranges with the same end point, keep the longest one (which is the first one with current sort order)
mymhends3() {
            i=$1;
            awk -F' ' -v i=$i '$5==i{print $0}' "$pd"/mhends2.tmp | sort -t' ' -k4,4n | head -1;
}
export -f mymhends3;

#myevalmicrohaps() takes the accounting of microhaplotype spans from mygetends() and calculates
#the number of alleles at each microhaplotype locus
myevalmicrohaps() {
                  i=$1; #incoming is "contig lesite resite lemhrange remhrange mhlength" like "51jcf7180007742276 11434 11434 11431 11440 10"
                 
                  contig=$(cut -d' ' -f1 <<<"$i"); #name of contig
                  col1=$(cut -d' ' -f6 <<<"$i"); #length of microhaploblock, may be used as number of columns for samtools tview to print
                  if [[ $col1 < 10 ]]; then col=10; else col="$col1"; fi; #shell variable $COLUMNS cannot be lower than 10
                  mhstart=$(cut -d' ' -f4 <<<"$i" | cut -d'-' -f1); #start position of microhaploblock 
                  mhend=$(cut -d' ' -f5 <<<"$i" | cut -d'-' -f1); #end position of microhaploblock 
                  
                  #use samtools tview to display microhaplotypes
                  mh=$(export COLUMNS="$col"; samtools tview -dT -p "$contig":"$mhstart" "$bam" | \
                      tail -n +4 | grep -v " " | awk -F' ' -v col1=$col1 '{print substr($1,1,col1)}' | \
                      sort | uniq -c | sed 's/^ *//');                  
                  #mh=$(export COLUMNS="$col"; /share/apps/samtools tview -dT -p "$contig":"$mhstart" "$bam" | \
                  #    tail -n +4 | grep -v " " | awk -F' ' -v col1=$col1 '{print substr($1,1,col1)}' | \
                  #    sort | uniq -c | sed 's/^ *//');                  
                  if [[ "$keepsingl" == NO ]];
                  then mh=$(grep -v ^"1 " <<<"$mh");
                  fi;
                  
                  #if no sequences are retained after removing singletons, or none exist, bounce out of this locus
                  #if [[ "$mh" == "" ]]; then return; fi;
                  if [[ "$mh" == "" ]];
                  then numseq=0; numalleles=0; mhcounts=0; freqs=0; alleleseqs="";
                  else mh=$(echo "$mh" | sort -t' ' -k1,1nr); #sort $mh to present alleles by decreasing frequency
                  
                    #assemble output string and report to parallel statement
                    numseq=$(cut -d' ' -f1 <<<"$mh" | awk '{s+=$1} END {print s}'); #total number of sequences considered
                  
                    numalleles=$(wc -l <<<"$mh"); #number of microhaplotype alleles
                    mhcounts=$(cut -d' ' -f1 <<<"$mh" | tr '\n' ':' | sed 's/:$//'); #counts of each allele
                    freqs=$(cut -d' ' -f1 <<<"$mh" | awk -v numseq=$numseq '{print $1/numseq}' | tr '\n' ':' | sed 's/:$//'); #frequencies of each allele
                    alleleseqs=$(cut -d' ' -f2 <<<"$mh" | tr '\n' ':' | sed 's/:$//'); #sequences of each allele
                  fi;
                  echo '#'"$contig $mhstart $mhend $col1 $numseq $numalleles $mhcounts $freqs $alleleseqs"; #report result to parallel statement
}
export -f myevalmicrohaps;

#mydd1() is used to rapidly collect full hapxmlog.txt output line for each mh with a unique start point, store it in $dd1
mydd1() {
       i=$1;
       awk -F' ' -v i=$i '$2==i{print $0}' "$pd"/a.tmp;
}
export -f mydd1;

#mydd2() is used to rapidly collect full hapxmlog.txt output line for each mh with a unique end point, store it in $dd2
mydd2() {
       i=$1;
       awk -F' ' -v i=$i '$3==i{print $0}' "$pd"/z.tmp;
}
export -f mydd2;

#mycc1(): for mhs with the same mhstart, find the mh with the most variants ($6 numalleles), if tie use max depth ($5, numseq),
#if tie use max length ($4, maxlength), if tie choose at random
mycc1() {
        i=$1;
        e=$(awk -F' ' -v i=$i '$2==i{print $0}' "$pd"/a.tmp | sort -t' ' -k6,6nr -k5,5nr -k4,4nr); #collect and sort full hapxmlog.txt output lines for current mhstart
      
        g=$(cut -d' ' -f4-6 <<<"$e" | head -1 | sed 's/^ *//'); #get fields 4-6 for the top entry
        j=$(cut -d' ' -f1 <<<"$g"); #value $4 mhlength
        k=$(cut -d' ' -f2 <<<"$g"); #value $5 numseq
        l=$(cut -d' ' -f3 <<<"$g"); #value $6 numalleles
        h=$(awk -F' ' -v i=$i -v j=$j -v k=$k -v l=$l '$2==i && $4==j && $5==k && $6==l{print $0}' <<<"$e"); #for current mhstart, collect lines that have identical fields 4-6 as top entry

        if (( $(wc -l <<<"$h") > 1 ));
        then shuf <<<"$h" | head -1; #if $h has more than one line choose one line randomly, echo it
        else echo "$h";
        fi;     
}
export -f mycc1;
    
#mycc2(): for mhs with the same mhend, find the mh with the most variants ($6 numalleles), if tie use max depth ($5, numseq),
#if tie use max length ($4, maxlength), if tie choose at random
mycc2() {
        i=$1;
        e=$(awk -F' ' -v i=$i '$3==i{print $0}' "$pd"/z.tmp | sort -t' ' -k6,6nr -k5,5nr -k4,4nr); #collect and sort full hapxmlog.txt output lines for current mhend
      
        g=$(cut -d' ' -f4-6 <<<"$e" | head -1 | sed 's/^ *//'); #get fields 4-6 for the top entry
        j=$(cut -d' ' -f1 <<<"$g"); #value $4 mhlength
        k=$(cut -d' ' -f2 <<<"$g"); #value $5 numseq
        l=$(cut -d' ' -f3 <<<"$g"); #value $6 numalleles
        h=$(awk -F' ' -v i=$i -v j=$j -v k=$k -v l=$l '$3==i && $4==j && $5==k && $6==l{print $0}' <<<"$e"); #for current mhend, collect lines that have identical fields 4-6 as top entry

        if (( $(wc -l <<<"$h") > 1 ));
        then shuf <<<"$h" | head -1; #if $h has more than one line choose one line randomly, echo it
        else echo "$h";
        fi;     
}
export -f mycc2;
  
#compute an input file for -u userrange if requested
mymhends4() {
            i=$1;
            k=$(cut -d: -f1 <<<"$i"); #mhstart
            l=$(cut -d: -f2 <<<"$i"); #mhend
            awk -F' ' -v k="$k" -v l="$l" '$4==k && $5==l{print $0}' "$pd"/mhends1.tmp; #echo the line corresponding to the the mh range retained in the variable tiling array
}
export -f mymhends4;

### END SUBROUTINES ###


#define variables and establish defaults
#bam, path to the bam file containing the alignment of reads to ref
#sites, genomic regions to use
stF=2048; #samtools view -F option, see https://broadinstitute.github.io/picard/explain-flags.html
stq=1; #samtools view -q option
debug=NO;
keepsingl=NO;
useuserranges=NO;
suppar=""; #suppress parallel contig extraction, default is allow GNU parallel --jobs equal to max, -sp switch will set $suppar to --jobs=1 for all parallel statements
tilearry=NO;
vartarry=NO;

#acquire command line variables
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b)
    bam="$2"
    shift # past argument
    shift # past value
    ;;
    -o)
    outfol="$2"
    shift # past argument
    shift # past value
    ;;
    -s)
    sites="$2"
    shift # past argument
    shift # past value
    ;;
    -F)
    stF="$2"
    shift # past argument
    shift # past value
    ;;
    -q)
    stq="$2"
    shift # past argument
    shift # past value
    ;;
    -ssh)
    ssh1="--sshloginfile $2"
    shift # past argument
    shift # past value
    ;;
    -u)
    useuserranges="YES"
    userrangefile="$2"
    shift # past argument
    shift # past value
    ;;
    -db)
    debug=YES
    shift # past argument
    ;;
    -ks)
    keepsingl=YES
    shift # past argument
    ;;
    -sp)
    suppar="--jobs 1";
    shift # past argument
    ;;
    -ta)
    tilearry=YES;
    shift # past argument
    ;;
    -va)
    vartarry=YES;
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

pd=$(pwd)"/$outfol"; export pd; #path to working directory
if [ ! -d "$pd" ]; then mkdir "$pd"; fi; #make the working directory if not already existing

e=$(cat "$sites"); #content of file $sites
if [[ "$useuserranges" == "YES" ]];
then u=$(cat "$userrangefile"); #content of file $userrange
fi;

export stF; 
export stq;
export bam;
export debug;
export keepsingl;

#log
log="$outfol"/hapxmlog.txt;
date > "$log";
echo >> "$log";
echo "Executable: $0" >> "$log";
echo "Alignment file (bam): $bam" >> "$log";
echo "Target sites: $sites" >> "$log";
echo "samtools view -F $stF" >> "$log";
echo "samtools view -q $stq" >> "$log";
echo "Use user-defined ranges: $useuserranges" >> "$log";
if [[ "$useuserranges" == "YES" ]]; then echo "User-defined ranges: $userrangefile" >> "$log"; fi;
echo "Keep singletons: $keepsingl" >> "$log";
echo "Debug on: $debug" >> "$log";
echo "Process short tiling array: $tilearry" >> "$log";
echo "Process variant-rich tiling array: $vartarry" >> "$log";
echo >> "$log";

#calculate ranges of microhaploblocks at contigname:site-range
mhends=$(echo "$e" | parallel $ssh1 $suppar --env stq --env stF --env bam --env debug --env mygetends mygetends);
if [[ "$debug" == "YES" ]]; then echo "$mhends" > "$pd"/mhends.txt; fi;
#sort on contig X microhaploblock range left end, then on unique microhaploblock ranges
mhends1=$(sed 's/[:-]/ /g' <<<"$mhends" | sort -t' ' -k1,1 -k4,4n | sort -u -t' ' -k1,1 -k4,4n -k5,5n);
if [[ "$debug" == "YES" ]]; then echo "$mhends1" > "$pd"/mhendssorted.txt; fi;
echo "$mhends1" > "$pd"/mhends1.tmp; #create a temp file for faster parallel access

#extract longest haploblocks across the tiling array
#for microhaploblock ranges with the same start point, keep the longest one (which is the last with current sort order)
#mhrstart=$(cut -d' ' -f4 <<<"$mhends1" | sort -un | tr '\n' ' '); #unique microhaploblock range start
#mhends2=$(for i in $mhrstart;
#  do awk -F' ' -v i=$i '$4==i{print $0}' <<<"$mhends1" | sort -t' ' -k5,5n | tail -1;
#  done;)
mhrstart=$(cut -d' ' -f4 <<<"$mhends1" | sort -un); #unique microhaploblock range start
mhends2=$(echo "$mhrstart" | parallel --env pd mymhends2);
mhends2=$(echo "$mhends2" | sort -t' ' -k1,1 -k4,4n | sort -u -t' ' -k1,1 -k4,4n -k5,5n);
echo "$mhends2" >  "$pd"/mhends2.tmp; #create a temp file for faster parallel access

#for microhaploblock ranges with the same end point, keep the longest one (which is the first one with current sort order)
#mhrend=$(cut -d' ' -f5 <<<"$mhends2" | sort -un | tr '\n' ' '); #unique microhaploblock range end
#mhends3=$(for i in $mhrend;
#  do awk -F' ' -v i=$i '$5==i{print $0}' <<<"$mhends2" | sort -t' ' -k4,4n | head -1;
#  done;)
#mhends3=$(sort -t' ' -k1,1 -k4,4n <<<"$mhends3" | sort -u -t' ' -k1,1 -k4,4n -k5,5n); #sort $mhends3 nicely
mhrend=$(cut -d' ' -f5 <<<"$mhends2" | sort -un); #unique microhaploblock range end
mhends3=$(echo "$mhrend" | parallel --env pd mymhends3);
mhends3=$(sort -t' ' -k1,1 -k4,4n <<<"$mhends3" | sort -u -t' ' -k1,1 -k4,4n -k5,5n); #sort $mhends3 nicely


#if user has supplied microhaploblock ranges by invoking the -u option, substitute those
#if user has elected to only process the short tiling array (-ta) use $mhends3
#otherwise process all unique mhs
if [[ "$useuserranges" == "YES" ]];
then mhendsin="$u";
elif [[ "$tilearry" == "YES" ]];
then mhendsin="$mhends3"; #short tiling array across all mhs
  if [[ "$debug" == "YES" ]]; then echo "$mhends3" > "$pd"/mhendstiled.txt; fi;
else mhendsin="$mhends1"; #all unique mh ranges, same as mhendssorted.txt
fi;

#count microhaploblock alleles at specified microhaploblock loci
echo "#contig mhstart mhend mhlength numseq numalleles counts freqs alleleseqs" >> "$log"; #add header for output table to log
report=$(echo "$mhendsin" | parallel $ssh1 $suppar --env bam --env debug \
    --env keepsingl --env myevalmicrohaps myevalmicrohaps);
report=$(sort -t' ' -k1,1 -k2,2n <<<"$report"); #sort by mh start position
echo "$report" >> "$log";



#if user has elected to calculate and process a variable tiling array (-va), extract those from hapxmlog.txt,
#identify mhs belonging to the variable tiling array and save to hapxmlogvar.txt
#save the variable tiling array in -u userrange format when $debug==YES 
if [[ "$vartarry" == "YES" ]];
then logv="$outfol"/hapxmlogvar.txt;
  grep -v ^'#' "$log" > "$logv"; #add common header from hapxmlog.txt
  echo "#contig mhstart mhend mhlength numseq numalleles counts freqs alleleseqs" >> "$logv"; #add header for output table to log

  a=$(grep ^'#' "$log" | awk -F' ' '$5>0{print $0}' | tail -n +2); #extract all lines with at least one mh called

  #filter to retain only one mh per mhstart point
  b=$(cut -d' ' -f2 <<<"$a" | sort -n | uniq -c | sed 's/^ *//'); #count the number of mhs per mhstart point
  c=$(grep -v ^"1 " <<<"$b" | cut -d' ' -f2);  #find mhstart with multiple mhs
  d=$(grep ^"1 " <<<"$b" | cut -d' ' -f2);  #find mhstart with only one mh
  
  #collect full hapxmlog.txt output line for each mh with a unique start point
  echo "$a" > "$pd"/a.tmp; #save a temporary copy of large variable $a
  dd1=$(echo "$d" | parallel --env pd mydd1);
  dd1=$(echo "$dd1" | sort -t' ' -k2,2n); #sort numeric on start point

  #for mhs with the same mhstart, find the mh with the most variants ($6 numalleles), if tie use max depth ($5, numseq),
  #if tie use max length ($4, maxlength), if tie choose at random
#  cc1=$(for i in $c;
#    do e=$(awk -F' ' -v i=$i '$2==i{print $0}' "$outfol"/a.tmp | sort -t' ' -k6,6nr -k5,5nr -k4,4nr); #collect and sort full hapxmlog.txt output lines for current mhstart
#      
#      g=$(cut -d' ' -f4-6 <<<"$e" | head -1 | sed 's/^ *//'); #get fields 4-6 for the top entry
#      j=$(cut -d' ' -f1 <<<"$g"); #value $4 mhlength
#      k=$(cut -d' ' -f2 <<<"$g"); #value $5 numseq
#      l=$(cut -d' ' -f3 <<<"$g"); #value $6 numalleles
#      h=$(awk -F' ' -v j=$j -v k=$k -v l=$l '$4==j && $5==k && $6==l{print $0}' "$outfol"/a.tmp); #collect full hapxmlog.txt output line(s) matching fields 4-6 for top entry
#
#      if (( $(wc -l <<<"$h") > 1 ));
#      then shuf <<<"$h" | head -1; #if $h has more than one line choose one line randomly, echo it
#      else echo "$h";
#      fi;     
#    done;)

  cc1=$(echo "$c" | parallel --env pd mycc1);

  z=$(echo "$cc1"$'\n'"$dd1" | sort -t' ' -k2,2n -k3,3n); #transfer lines filtered for mhstart to starting variable
  echo "$z" > "$pd"/z.tmp; #save a temporary copy of large variable $z

  #starting with mhs with unique mhstart points, filter to retain only one mh per mhend point
  b2=$(cut -d' ' -f3 <<<"$z" | sort -n | uniq -c | sed 's/^ *//'); #count the number of mhs per mhend point
  c2=$(grep -v ^"1 " <<<"$b2" | cut -d' ' -f2);  #find mhend with multiple mhs
  d2=$(grep ^"1 " <<<"$b2" | cut -d' ' -f2);  #find mhend with only one mh
  
  #collect full hapxmlog.txt output line for each mh with a unique end point
  dd2=$(echo "$d2" | parallel --env pd mydd2);
  dd2=$(echo "$dd2" | sort -t' ' -k2,2n); #sort numeric on start point
  
#  dd2=$(for i in $d;
#         do awk -F' ' -v i=$i '$3==i{print $0}' <<<"$z";
#         done;)
         
  #for mhs with the same mhend, find the mh with the most variants ($6 numalleles), if tie use max depth ($5, numseq),
  #if tie use max length ($4, maxlength), if tie choose at random
#  cc2=$(for i in $c2;
#    do e=$(awk -F' ' -v i=$i '$3==i{print $0}' "$outfol"/z.tmp | sort -t' ' -k6,6nr -k5,5nr -k4,4nr); #collect and sort full hapxmlog.txt output lines for current mhend
#      
#      g=$(cut -d' ' -f4-6 <<<"$e" | head -1 | sed 's/^ *//'); #get fields 4-6 for the top entry
#      j=$(cut -d' ' -f1 <<<"$g"); #value $4 mhlength
#      k=$(cut -d' ' -f2 <<<"$g"); #value $5 numseq
#      l=$(cut -d' ' -f3 <<<"$g"); #value $6 numalleles
#      h=$(awk -F' ' -v j=$j -v k=$k -v l=$l '$4==j && $5==k && $6==l{print $0}' "$outfol"/z.tmp); #collect full hapxmlog.txt output line(s) matching fields 4-6 for top entry
#
#      if (( $(wc -l <<<"$h") > 1 ));
#      then shuf <<<"$h" | head -1; #if $h has more than one line choose one line randomly, echo it
#      else echo "$h";
#      fi;     
#    done;)

  cc2=$(echo "$c2" | parallel --env pd mycc2);

  #report to log
  zz=$(echo "$cc2"$'\n'"$dd2" | sort -t' ' -k2,2n -k3,3n | awk 'NF');
  echo "$zz" >> "$logv";

  #compute an input file for -u userrange if requested
  if [[ "$debug" == "YES" ]];
  then j=$(cut -d' ' -f2-3 <<<"$zz" | tr ' ' ':');#get list of mh ranges for variable tiling array
    mhends4=$(echo "$j" | parallel --env pd mymhends4);
    
    #    mhends4=$(for i in $j;
#                do k=$(cut -d: -f1 <<<"$i"); #mhstart
#                  l=$(cut -d: -f2 <<<"$i"); #mhend
#                  awk -F' ' -v k="$k" -v l="$l" '$4==k && $5==l{print $0}' <<<"$mhends1"; #echo the line corresponding to the the mh range retained in the variable tiling array
#                done;)
    
    echo "$mhends4" > "$pd"/mhendsvar.txt;
  fi;
  
  #clean up
    rm "$pd"/mhends1.tmp;
    rm "$pd"/mhends2.tmp;
    rm "$pd"/a.tmp; #remove temporary copy of large variable $a
    rm "$pd"/z.tmp; #remove temporary copy of large variable $z

fi; #vartarry






