READ1 = /ebio/abt3_projects2/hic_metagenomics_pilot/mbermejo/TEST/ecoli_microc/WT_netropsin_20/SRR33290678_1.fastq.gz
READ2 = /ebio/abt3_projects2/hic_metagenomics_pilot/mbermejo/TEST/ecoli_microc/WT_netropsin_20/SRR33290678_2.fastq.gz
FASTA = /ebio/abt3_projects2/hic_metagenomics_pilot/mbermejo/TEST/ecoli_microc/MG1655.fasta
BS = 1000
OUTDIR = /tmp/global2/mbermejo/conda/hicexplorer/Micro-c

#get allignments of the microC reads against the reference genomes
bwa mem -A1 -B4 -E50 -L0 -t 20 $FASTA $READ1 | samtools view -@ 20 -Shb - > $OUTDIR/mate_R1.bam
bwa mem -A1 -B4 -E50 -L0 -t 20 $FASTA $READ2 | samtools view -@ 20 -Shb - > $OUTDIR/mate_R2.bam

#Calculate contact matrix for a specific bin size
hicBuildMatrixMicroC --samFiles $OUTDIR/mate_R1.bam $OUTDIR/mate_R2.bam --outFileName $OUTDIR/${BS}b/contact_matrix.h5 --QCfolder $OUTDIR/${BS}b/QC --threads 8 --genomeAssembly $FASTA -bs $BS
#Normalize the matrix
hicCorrectMatrix diagnostic_plot -m $OUTDIR/${BS}b/contacts_matrix.h5 -o $OUTDIR/${BS}b/contacts_plot.png
hicCorrectMatrix correct -m $OUTDIR/${BS}b/contacts_matrix.h5 --filterThreshold -1 5 -o $OUTDIR/${BS}b/contacts_matrix_corrected.h5

#Change matrix format
hicConvertFormat --matrices $OUTDIR/${BS}b/contacts_matrix_corrected.h5 --outFileName $OUTDIR/${BS}b/contact_matrix --inputFormat h5 --outputFormat ginteractions
hicConvertFormat -m 5000b/contact_matrix_corrected.h5 1000b/contact_matrix_corrected.h5 900b/contact_matrix_corrected.h5 800b/contact_matrix_corrected.h5 700b/contact_matrix_corrected.h5 600b/contact_matrix_corrected.h5 500b/contact_matrix_corrected.h5 400b/contact_matrix_corrected.h5 300b/contact_matrix_corrected.h5 200b/contact_matrix_corrected.h5 100b/contact_matrix_corrected.h5 50b/contact_matrix_corrected.h5 20b/contact_matrix_corrected.h5 --inputFormat h5 --outputFormat mcool -o all_resolutions.mcool # --resolutions 5000 1000 900 800 700 600 500 400 300 200 100 50 20 

#Plot heatmaps
hicPlotMatrix -m $OUTDIR/${BS}b/contact_matrix_corrected.h5 -out $OUTDIR/${BS}b/heatmap_full.png --log --dpi 300
hicPlotMatrix -m $OUTDIR/${BS}b/contact_matrix_corrected.h5 -out $OUTDIR/${BS}b/heatmap-2350-2600.png --log --region NC_000913.3:2350000-2600000 --colorMap hot_r --dpi 200
hicPlotMatrix -m $OUTDIR/${BS}b/contact_matrix_corrected.h5 -out $OUTDIR/${BS}b/heatmap_E_faecium.png --log --chromosomeOrder CP118955.1 CP118956.1  --dpi 300