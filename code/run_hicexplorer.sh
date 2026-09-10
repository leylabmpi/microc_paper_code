READ1 = R1.fastq.gz #Micro-C R1 file
READ2 = R2.fastq.gz #Micro-C R2 file
FASTA = reference.fasta #Fasta file containing reference genomes
BS = 1000 #bin size
OUTDIR = Hicexplorer_out #output directory name

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

#Plot heatmaps
hicPlotMatrix -m $OUTDIR/${BS}b/contact_matrix_corrected.h5 -out $OUTDIR/${BS}b/heatmap_full.png --log --dpi 300
hicPlotMatrix -m $OUTDIR/${BS}b/contact_matrix_corrected.h5 -out $OUTDIR/${BS}b/heatmap-2350-2600.png --log --region NC_000913.3:2350000-2600000 --colorMap hot_r --dpi 200
hicPlotMatrix -m $OUTDIR/${BS}b/contact_matrix_corrected.h5 -out $OUTDIR/${BS}b/heatmap_E_faecium.png --log --chromosomeOrder CP118955.1 CP118956.1  --dpi 300