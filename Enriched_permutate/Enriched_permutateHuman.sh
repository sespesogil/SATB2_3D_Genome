
in=$1
roi1=$2
roi2=$3
res=$4
it=$5

echo "Hello stranger... welcome"

if [[ -f $in/ ]]
then
    echo "Welcome back stranger ;-)"
   # rm -r $in/consensus $in/toUPset $in/intersection_annotations $in/aggregate.ann.R $in/toUPset.R 
fi

module load gcc/9.3.0-gcc-9.1.0-dbq7v4m

case $4 in
        0.5kb )
        univ=/resources/human/genome.human.0.5kb.bed
        ;;
        500bp)
        univ=/resources/human/genome.human.0.5kb.bed
        ;;
        1kb )
        univ=/resources/human/genome.human.1kb.bed
        ;;
        5kb )
        univ=/resources/human/genome.human.5kb.bed
        ;;
        10kb )
        univ=/resources/human/genome.human.10kb.bed
        ;;
        25kb )
        univ=/resources/human/genome.human.25kb.bed
        ;;
        50kb )
        univ=/human/genome.human.50kb.bed
        ;;
        100kb )
        univ=/genome.human.100kb.bed
        ;;
        * )
        echo "Wrong resolution provided. Available resolutions: 0.5kb, 1kb, 5kb, 10kb, 25kb, 50kb, 100kb" 
        esac

# permutation:

cat > $in/permutation.R << EOF
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dir<-args[1]
roi1 <- args[2]
roi2<- args[3]
univ <- args[4]
iterations = as.numeric(args[5])

require(regioneR)

targets <- toGRanges(paste(dir,roi1, sep="/"))
regions <- toGRanges(paste(dir,roi2, sep="/"))
consensus<- toGRanges(univ)

pt <- permTest(A=targets, ntimes=iterations, randomize.function=resampleRegions, universe=consensus, evaluate.function=numOverlaps, B=regions, verbose=FALSE)
summary(pt)

sink(paste0(dir, "/summary.txt"))
summary(pt)
sink()

pdf(paste0("permutation.pdf"))
  plot(pt)
  dev.off()
EOF

# launcher

echo "Launching permutation analysis... "


cat > $in/permutationTOslurm.slrm << EOF
#!/bin/bash
##
#SBATCH -J Enriched_perm
#SBATCH -n 6
#SBATCH --mem=60G
#SBATCH --ntasks-per-node=40

in_dir=$1
roi1=$2
roi2=$3
univ=$4
it=$5

module load gcc/9.1.0-gcc-4.8.5-mj7s6dg

Rscript $1/permutation.R $1 $roi1 $roi2 $univ $it

EOF


sbatch $in/permutationTOslurm.slrm $in $roi1 $roi2 $univ $it

echo "...done
