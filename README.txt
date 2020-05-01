hapxm identifies microhaploblock variation in bam files

In working folder:
1) A bam file containing sequences mapped by 'hapx -mb'
2) A line delimited list of target sites

Requirements (in path):
1) samtools (hapxm calls samtools view/tview)
3) GNU parallel

Usage: hapxm -b bam -o out [-F exc] [-q qual] [-db -ks] -s sites
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

-db = debugging mode, save some internal data structures as files (may produce a lot of files)
-ks = keep singletons, default behavior is to ignore microhaplotype singletons (occur in only 1 sequence)

Examples: hapxm.sh -b /share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/Hs1pro1l1.finalaln.bam \
            -o hxm1 -db -s <(for i in $(seq 11423 1 11500); do echo 51jcf7180007742276:"$i"-"$i"; done;)

hapxm.sh -b /share/space/reevesp/patellifolia/hapxtest/hapxsummary/bwamem/Hs1pro1l1.finalaln.bam \
            -o hxm1 -db -s <(for i in $(seq 9889 1 15245); do echo 51jcf7180007742276:"$i"-"$i"; done;)

