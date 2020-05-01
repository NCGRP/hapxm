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
              g=$(/share/apps/samtools view -F "$stF" -q "$stq" "$bam" "$i" | cut -d$'\t' -f4,6); #-F 2048 excludes supplementary alignments
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
              
              
              
#              YOU ARE HERE, working on expanding the cigar string
#              echo $cig | sed 's:\([DIMSX=]\):\1\n:g' | sed '/^$/d' | awk -F[[:alpha:]] '{print $1}'
}
export -f mygetends;

#myevalmicrohaps() takes the accounting of microhaplotype spans from mygetends() and calculates
#the number of alleles at each microhaplotype locus
myevalmicrohaps() {
                  i=$1; #incoming is "contig:site-range microhaploblockrange microhaploblocklength" like "51jcf7180007742276:11500-11500 11498-11500 3"
                 
                  contig=$(cut -d: -f1 <<<"$i"); #name of contig
                  col1=$(cut -d' ' -f3 <<<"$i"); #length of microhaploblock, may be used as number of columns for samtools tview to print
                  if [[ $col1 < 10 ]]; then col=10; else col=col1; fi; #shell variable $COLUMNS cannot be lower than 10
                  mhstart=$(cut -d' ' -f2 <<<"$i" | cut -d'-' -f1); #start position of microhaploblock 
                  
                  
                  
                  
                  #export COLUMNS="$col"; /share/apps/samtools tview -dT -p "$contig":"$mhstart" "$bam" | tail -n +4 | grep -v " " | sort | uniq -c;                  





                  #export COLUMNS="$col"; /share/apps/samtools tview -dC -p "$contig":"$mhstart" "$bam" | tail -n +4 | grep -v " " | sort | uniq -c;                  
}
export -f myevalmicrohaps;

### END SUBROUTINES ###




#define variables and establish defaults
#bam, path to the bam file containing the alignment of reads to ref
#sites, genomic regions to use
stF=2048; #samtools view -F option, see https://broadinstitute.github.io/picard/explain-flags.html
stq=1; #samtools view -q option

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
    -db)
    debug=YES
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
export stF; 
export stq;
export bam;

#log
log="$outfol"/hapxmlog.txt;
date > "$log";
echo >> "$log";
echo "Executable: $0" >> "$log";
echo "Alignment file (bam): $bam" >> "$log";
echo "Target sites: $sites" >> "$log";
echo "samtools view -F $stF" >> "$log";
echo "samtools view -q $stq" >> "$log";
echo >> "$log";

#calculate ranges of microhaploblocks at contigname:site-range
mhends=$(echo "$e" | parallel --sshloginfile /home/reevesp/machines --env stq --env stF --env bam --env mygetends mygetends);
if [[ "$debug" == "YES" ]]; then echo "$mhends" > "$pd"/mhends.txt; fi;
#sort on contig X microhaploblock range left end, then on unique microhaploblock ranges
mhends1=$(sed 's/[:-]/ /g' <<<"$mhends" | sort -t' ' -k1,1 -k4,4n | sort -u -t' ' -k4,5);

#extract longest haploblocks across the tiling array
#for microhaploblock ranges with the same start point, keep the longest one (which is the last with current sort order)
mhrstart=$(cut -d' ' -f4 <<<"$mhends1" | sort -u | tr '\n' ' '); #unique microhaploblock range start
mhends2=$(for i in $mhrstart;
  do awk -F' ' -v i=$i '$4==i{print $0}' <<<"$mhends1" | tail -1;
  done;)

#for microhaploblock ranges with the same start point, keep the longest one (which is the first one with current sort order)
mhrend=$(cut -d' ' -f5 <<<"$mhends2" | sort -u | tr '\n' ' '); #unique microhaploblock range end
mhends3=$(for i in $mhrend;
  do awk -F' ' -v i=$i '$5==i{print $0}' <<<"$mhends2" | head -1;
  done;)


mhends=$(sort -t':' -k1,1 <<<"$mhends" | sort -t'-' -k2,2n | sort -u -t' ' -k2,2 | awk -F[' '-:] '{print $1, $4, $0}'| sort -t' ' -k1,1 -k2,2n | cut -d' ' -f3-);
mhends=$(echo "$mhends" | sed 's/[:-]/ /g' | sort -t' ' -u -k1,1 -k4,5n -uk4,5 -k4,4n |           sort -t':' -k1,1 <<<"$mhends" | sort -t'-' -k2,2n | sort -u -t' ' -k2,2 | awk -F[' '-:] '{print $1, $4, $0}'| sort -t' ' -k1,1 -k2,2n | cut -d' ' -f3-);

echo "$mhends" | sed 's/[:-]/ /g' | sort -t' ' -k1,1 -k2,2n | sort -u -t' ' -k4,5 | md5sum
echo "$mhends" | sed 's/[:-]/ /g' | sort -t' ' -k1,1 -k4,4n | sort -u -t' ' -k4,5 | md5sum

if [[ "$debug" == "YES" ]]; then echo "$mhends" > "$pd"/mhendssorted.txt; fi;
#determine path through contig that maximizes


#count microhaploblock alleles at each microhaploblock locus
echo "$mhends" | parallel --sshloginfile /home/reevesp/machines --env bam --env myevalmicrohaps myevalmicrohaps;



#    YOU ARE HERE
#    #called like: 
#    #mhends=$( (for i in $(seq 11423 1 11500); do echo 51jcf7180007742276:"$i"-"$i"; done;) \
#    #  | parallel --sshloginfile /home/reevesp/machines --env stq --env bam --env mygetends mygetends );
#    #mhends=$(sort -t'-' -k1,1n <<<"$mhends");
#    
#    
#    
#    
#    pd=$(pwd);
#    bam=/share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/Hs1pro1l1.finalaln.bam;
#    export bam;
#    
#    mhends=$( (for i in $(seq 11423 1 11500); do echo 51jcf7180007742276:"$i"-"$i"; done;) \
#      | parallel --sshloginfile /home/reevesp/machines --env stq --env bam --env mygetends mygetends );
#    mhends=$(sort -t'-' -k1,1n <<<"$mhends");
#    
#    echo "$mhends" | parallel --sshloginfile /home/reevesp/machines --env bam --env myevalmicrohaps myevalmicrohaps;
#    
#    
#    
#    export COLUMNS=200; samtools tview -dT -p 51jcf7180007742276:11491 "$bam" | tail -n +4 | grep -v " " | sort | uniq -c;
#    
#    
#    /share/apps/samtools view -F 2048 -q "$stq" -O BAM "$bam" "$i" 2>/dev/null > "$i".TMP.bam;
#    samtools index "$i".TMP.bam 2>/dev/null;
#    export COLUMNS=10; samtools tview -dT -p 51jcf7180007742276:11434 "$i".TMP.bam | tail -n +4;
#    
#    export COLUMNS=10; samtools tview -dT -p 51jcf7180007742276:11434 "$i".TMP.bam | tail -n +4 | sort | uniq -c;
#    export COLUMNS=10; samtools tview -dT -p 51jcf7180007742276:11431 "$i".TMP.bam | tail -n +4 | sort | uniq -c;
#    export COLUMNS=10; samtools tview -dT -p 51jcf7180007742276:11434 "$bam" | tail -n +4 | sort | uniq -c;
#    export COLUMNS=10; samtools tview -dT -p 51jcf7180007742276:11431 "$bam" | tail -n +4 | sort | uniq -c;
#    
#    export COLUMNS=10; samtools tview -dT -p 51jcf7180007742276:11431 "$bam" | tail -n +4 | grep -v " " | sort | uniq -c;
#                   ^length from $mhends                         ^le from $mhends
#    export COLUMNS=200; samtools tview -dT -p 51jcf7180007742276:11491 "$bam" | tail -n +4 | grep -v " " | sort | uniq -c;
#    
#    
#    i=51jcf7180007742276:11434-11434;
#    ii=51jcf7180007742276:11433-11433;
#    /share/apps/samtools view -F 2048 -q "$stq" "$bam" "$i" 2>/dev/null > "$i".TMP.sam;
#    /share/apps/samtools view -F 2048 -q "$stq" "$bam" "$ii" 2>/dev/null > "$ii".TMP.sam;
#    
#    
#    /share/apps/samtools view -q "$stq" -O BAM "$bam" "$i" 2>/dev/null | \
#      /share/apps/samtools fasta - 2>/dev/null | \
#      /share/apps/bwa mem "$ref" - 2>/dev/null | \
#      /share/apps/samtools sort -O BAM 2>/dev/null | \
#      /share/apps/samtools view -F 2048 -q "$stq" -O BAM - "$i" 2>/dev/null | 
#     This ^ not working. Trying to extract reads for le-re region, realign, then print microhaplotype with tview
#    
#    
#    
#    
#    samtools tview -dT -p 51jcf7180007742276:11434 "$bam" | tail -n +3;
#    samtools tview -dC -p 51jcf7180007742276:11434 "$i".TMP.bam | tail -n +4;
#    
#    export COLUMNS=10; samtools tview -dT --reference "$ref" -p 51jcf7180007742276:11434 "$bam" | less
#    export COLUMNS=10; samtools tview -dT -p 51jcf7180007742276:11434 "$bam" | tail -n +3 | grep -v ^" " 
#    
 
