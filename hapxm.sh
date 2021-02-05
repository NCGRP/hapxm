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
echo >> "$log";

#calculate ranges of microhaploblocks at contigname:site-range
mhends=$(echo "$e" | parallel $ssh1 $suppar --env stq --env stF --env bam --env debug --env mygetends mygetends);
if [[ "$debug" == "YES" ]]; then echo "$mhends" > "$pd"/mhends.txt; fi;
#sort on contig X microhaploblock range left end, then on unique microhaploblock ranges
mhends1=$(sed 's/[:-]/ /g' <<<"$mhends" | sort -t' ' -k1,1 -k4,4n | sort -u -t' ' -k1,1 -k4,4n -k5,5n);
if [[ "$debug" == "YES" ]]; then echo "$mhends1" > "$pd"/mhendssorted.txt; fi;

#extract longest haploblocks across the tiling array
#for microhaploblock ranges with the same start point, keep the longest one (which is the last with current sort order)
mhrstart=$(cut -d' ' -f4 <<<"$mhends1" | sort -un | tr '\n' ' '); #unique microhaploblock range start
mhends2=$(for i in $mhrstart;
  do awk -F' ' -v i=$i '$4==i{print $0}' <<<"$mhends1" | sort -t' ' -k5,5n | tail -1;
  done;)

#for microhaploblock ranges with the same end point, keep the longest one (which is the first one with current sort order)
mhrend=$(cut -d' ' -f5 <<<"$mhends2" | sort -un | tr '\n' ' '); #unique microhaploblock range end
mhends3=$(for i in $mhrend;
  do awk -F' ' -v i=$i '$5==i{print $0}' <<<"$mhends2" | sort -t' ' -k4,4n | head -1;
  done;)

if [[ "$debug" == "YES" ]]; then echo "$mhends3" | sort -t' ' -k1,1 -k4,4n | sort -u -t' ' -k1,1 -k4,4n -k5,5n > "$pd"/mhendstiled.txt; fi;

#if user has supplied microhaploblock ranges by invoking the -u option, substitute those
#for $mhends3 here
if [[ "$useuserranges" == "YES" ]];
then mhends3="$u";
fi;

#count microhaploblock alleles at minimal tiling path microhaploblock loci
echo "#contig mhstart mhend mhlength numseq numalleles counts freqs alleleseqs" >> "$log"; #add header for output table to log
report=$(echo "$mhends3" | parallel $ssh1 $suppar --env bam --env debug \
    --env keepsingl --env myevalmicrohaps myevalmicrohaps);
report=$(sort -t' ' -k1,1 -k2,2n <<<"$report"); #sort by mh start position
echo "$report" >> "$log";
