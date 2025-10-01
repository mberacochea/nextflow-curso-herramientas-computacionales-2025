bwa index reference.fasta

bwa mem reference.fasta test_1.fastq.gz test_1.fastq.gz > aligned.sam

samtools view -bS aligned.sam > aligned.bam

samtools sort aligned.bam -o aligned.sorted.bam

samtools index aligned.sorted.bam