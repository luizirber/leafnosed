import os
from snakemake.utils import makedirs, listfiles

include: "conf.py"
workdir: WORKDIR

localrules: download_raw_data

rule download_raw_data:
    output:
        ['inputs/' + i for i in RAW_DATA]
    shell:
	    "wget -nH -m -P inputs/ --cut-dirs=2 ftp://crocgenomes.org/pub/phyllostomid_genomes/"

# http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm
# <sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz 
 
rule fastqc_illumina:
    input: 'inputs/HiSeq/Run_1305/Project_hanson/Sample_{name}/{name}_{barcode}_L{lane,\d+}_R{read_number}_001.fastq.gz'
    output: 'outputs/fastqc_illumina/Sample_{name}/{name}_{barcode}_L{lane,\d+}_R{read_number}_fastqc.html'
    shell: 'fastqc --casava -o $(dirname {output}) {input}'


#rule fastqc_pacbio:
#    input: 'inputs/PacBio/{sample_number}/_A01_1_Cell1_PacBioRun167_JHanson560_TK163824/Analysis_Results/m150428_212227_42131_c100799742550000001823174309091570_s1_p0.1.subreads.fastq'
    #input: 'inputs/PacBio/TK163824/_A01_1_Cell1_PacBioRun167_JHanson560_TK163824/Analysis_Results/m150428_212227_42131_c100799742550000001823174309091570_s1_p0.1.subreads.fastq'


rule publish_fastqc:
    input:
        "outputs/fastqc_illumina/Sample_46394/46394_GCCAAT_L001_R1_fastqc.html",
        "outputs/fastqc_illumina/Sample_46394/46394_GCCAAT_L001_R2_fastqc.html",
        "outputs/fastqc_illumina/Sample_46394/46394_GCCAAT_L002_R1_fastqc.html",
        "outputs/fastqc_illumina/Sample_46394/46394_GCCAAT_L002_R2_fastqc.html",
        "outputs/fastqc_illumina/Sample_46395/46395_CTTGTA_L001_R1_fastqc.html",
        "outputs/fastqc_illumina/Sample_46395/46395_CTTGTA_L001_R2_fastqc.html",
        "outputs/fastqc_illumina/Sample_46395/46395_CTTGTA_L002_R1_fastqc.html",
        "outputs/fastqc_illumina/Sample_46395/46395_CTTGTA_L002_R2_fastqc.html",
    shell: """
        ssh {REMOTE_HOST} "mkdir -p {REMOTE_PATH}/bat"
        scp {input} {REMOTE_HOST}:{REMOTE_PATH}/bat"
    """
