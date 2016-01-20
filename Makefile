snakemake: inputs files
	snakemake publish_fastqc

inputs:
	wget -nH -m -P inputs/ --cut-dirs=2 ftp://crocgenomes.org/pub/phyllostomid_genomes/

files:
	cd inputs && \
	find -type f | cut -c3- | grep -v .listing | sort > ../files

#outputs/pacbio_lordec/desRot.fastq.gz
#outputs/pacbio_lordec/macCal.fastq.gz

KSIZE=21

dists: outputs/abundance_dist/corr.macCal.dist outputs/abundance_dist/corr.desRot.dist \
       outputs/abundance_dist/illumina.macCal.dist outputs/abundance_dist/illumina.desRot.dist
#outputs/abundance_dist/%.unique: outputs/pacbio_lordec/%

outputs/abundance_dist/%.unique: outputs/sga/%/illumina_merged.fastq.gz
	env OMP_NUM_THREADS=16 unique-kmers.py -R $@ -e 0.01 -k ${KSIZE} $^

#outputs/abundance_dist/%.kh: outputs/pacbio_lordec/% outputs/abundance_dist/%.unique
outputs/abundance_dist/%.kh: outputs/sga/%/illumina_merged.fastq.gz outputs/abundance_dist/%.unique
	load-into-counting.py -k ${KSIZE} $(shell python estimate.py $(shell cat $(subst kh,unique,$@) | head -1| cut -d " " -f1) 0.05) $@ $<

outputs/abundance_dist/corr.%.dist: outputs/abundance_dist/%.kh outputs/pacbio_lordec/%
	abundance-dist.py -s $^ $@

outputs/abundance_dist/illumina.%.dist: outputs/abundance_dist/%.kh outputs/sga/%/illumina_merged.fastq.gz
	abundance-dist.py -s $^ $@


.SECONDARY:
