module load python/3-Anaconda
source activate SWI_SNF
cd /dartfs-hpc/rc/home/9/f0029s9/Homer_analysis/Motif

perl /dartfs-hpc/rc/home/9/f0029s9/.conda/envs/SWI_SNF/share/homer-4.10-0/.//configureHomer.pl -install hg38

scanMotifGenomeWide.pl custom.motifs hg38 -bed > /dartfs-hpc/rc/home/9/f0029s9/Homer_analysis/Motif_scan/hg38_motif.scan