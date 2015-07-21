snakemake: inputs files
	snakemake publish_fastqc

inputs:
	wget -nH -m -P inputs/ --cut-dirs=2 ftp://crocgenomes.org/pub/phyllostomid_genomes/

files:
	cd inputs && \
	find -type f | cut -c3- | grep -v .listing | sort > ../files

.SECONDARY:
