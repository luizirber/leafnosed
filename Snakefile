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
    input:
        lambda w: 'inputs/HiSeq/Run_1305/Project_hanson/Sample_' + SPECIES_TO_SAMPLE[w.species] + '/{name}_{barcode}_L{lane,\d+}_R{read_number}_001.fastq.gz'
    output: 'outputs/fastqc_illumina/{species}/{name}_{barcode}_L{lane,\d+}_R{read_number}_fastqc.html'
    shell: 'fastqc --casava -o $(dirname {output}) {input}'

def sga_inputs(w):
    PATH = 'inputs/HiSeq/Run_1305/Project_hanson/'
    final = []
    inputs = []

    for f in ILLUMINA_RAW[w.species]:
        seq_id = 'Sample_{sample}/{name}_{barcode}_L{lane}_R1'.format(
            name=w.name,
            barcode=w.barcode,
            lane=w.lane,
            sample=SPECIES_TO_SAMPLE[w.species])
        if seq_id in f:
            inputs.append(f)

    for i in sorted(inputs):
        final.append(i)
        final.append(i.replace("R1", "R2"))

    return final

rule sga_preprocess_illumina:
    input: sga_inputs,
    output: 'outputs/sga/{species}/{name}_{barcode}_L{lane,\d+}.fastq.gz'
    shell: 'sga preprocess --pe-mode 1 {input} | gzip -c > {output}'

rule sga_index_illumina:
#    input: 'outputs/sga/{species}/{name}_{barcode}_L{lane,\d+}.fastq.gz'
    input: '{prefix}.fastq.gz'
    output:
        '{prefix}.bwt',
        '{prefix}.sai'
    threads: 16
    shell: '''
       cd $(dirname {input})
       sga index -a ropebwt --no-reverse -t {threads} $(basename {input})
    '''

rule sga_preprocess_pacbio:
    input: 'outputs/pacbio_filtered/{species}_filtered_subreads.fastq.gz'
    output: 'outputs/sga/{species}/filtered_subreads.fastq.gz'
    shell: 'sga preprocess --pe-mode 0 {input} | gzip -c > {output}'

rule sga_index_pacbio:
#    input: 'outputs/sga/{species}/{name}_{barcode}_L{lane,\d+}.fastq.gz'
    input: 'outputs/sga/{species}/filtered_subreads.fastq.gz'
    output:
        'outputs/sga/{species}/filtered_subreads.bwt',
        'outputs/sga/{species}/filtered_subreads.sai',
    threads: 16
    shell: '''
       cd $(dirname {input})
       sga index -a sais -d 400000 --no-reverse -t {threads} $(basename {input})
    '''

rule sga_preqc:
    input:
        reads='{prefix}.fastq.gz',
        bwt='{prefix}.bwt',
        sai='{prefix}.sai',
    output: '{prefix}.preqc'
    threads: 16
    shell: '''
        cd $(dirname {input.reads})
        sga preqc -t {threads} $(basename {input.reads}) > $(basename {output})
    '''

rule sga_preqc_report:
    input: '{prefix}.preqc'
    output: '{prefix}.pdf'
    shell: 'scripts/sga-preqc-report.py -o {wildcards.prefix} {input}'

#############################################################################

rule pacbio_filtering:
    input:
        params='inputs/pacbio_filtering/{species}/params.xml',
        inputs='outputs/pacbio_filtering/{species}/input.xml'
#    output: ''
    shell:
        'smrtpipe -D NPROC=3 -D CLUSTER=BASH -D MAX_THREADS=4 '
        '--params={WORKDIR}/{input.params} xml:{WORKDIR}/{input.inputs} '
        '--output {WORKDIR}/outputs/pacbio_filtering/'

rule pacbio_filtering_convert_input:
    input: 'outputs/pacbio_filtering/{species}/input.fofn'
    output: 'outputs/pacbio_filtering/{species}/input.xml'
    shell:
        'fofnToSmrtpipeInput.py {input} > {output}'

rule pacbio_filtering_fofn:
    input:
         lambda w: ['inputs/' + f for f in PACBIO_RAW[w.species]
                                if f.endswith(".bax.h5")]
    output: 'outputs/pacbio_filtering/{species}/input.fofn'
    run:
        with open(output[0], 'w') as fofn:
            fofn.write("\n".join(os.path.join(WORKDIR, f) for f in input))
            fofn.write("\n")

#############################################################################

rule fastqc_pacbio:
    input: 'outputs/pacbio_filtered/{species}_filtered_subreads.fastq.gz'
    output: 'outputs/fastqc_pacbio/{species}/{species}_filtered_subreads_fastqc.html'
    shell: 'fastqc -o $(dirname {output}) {input}'

#rule fastqc_pacbio:
#    input: 'inputs/PacBio/{sample_number}/_A01_1_Cell1_PacBioRun167_JHanson560_TK163824/Analysis_Results/m150428_212227_42131_c100799742550000001823174309091570_s1_p0.1.subreads.fastq'
    #input: 'inputs/PacBio/TK163824/_A01_1_Cell1_PacBioRun167_JHanson560_TK163824/Analysis_Results/m150428_212227_42131_c100799742550000001823174309091570_s1_p0.1.subreads.fastq'

rule lordec_merge_short:
    input:
        lambda w: [i.format(species=w.species,
                            name=SPECIES_TO_SAMPLE[w.species],
                            barcode=SPECIES_TO_BARCODE[w.species]) for i in
        ('outputs/sga/{species}/{name}_{barcode}_L001.fastq.gz',
        'outputs/sga/{species}/{name}_{barcode}_L002.fastq.gz')]
    output:
        'outputs/sga/{species}/illumina_merged.fastq.gz'
    shell:
        'cat {input} > {output}'

rule lordec_correct:
    input:
      long='outputs/pacbio_filtered/{species}_filtered_subreads.fastq.gz',
      short='outputs/sga/{species}/illumina_merged.fastq.gz'
    output: 'outputs/pacbio_lordec/{species}.fastq.gz'
    threads: 16
    shell:
        "lordec-correct -t 5 -b 200 -e 0.4 -T {threads} -k 15 -a 20 -s 3 -i {input.long} -2 {input.short} -o {output}"


rule publish_fastqc:
    input:
        "outputs/fastqc_illumina/macCal/46394_GCCAAT_L001_R1_fastqc.html",
        "outputs/fastqc_illumina/macCal/46394_GCCAAT_L001_R2_fastqc.html",
        "outputs/fastqc_illumina/macCal/46394_GCCAAT_L002_R1_fastqc.html",
        "outputs/fastqc_illumina/macCal/46394_GCCAAT_L002_R2_fastqc.html",
        "outputs/fastqc_illumina/desRot/46395_CTTGTA_L001_R1_fastqc.html",
        "outputs/fastqc_illumina/desRot/46395_CTTGTA_L001_R2_fastqc.html",
        "outputs/fastqc_illumina/desRot/46395_CTTGTA_L002_R1_fastqc.html",
        "outputs/fastqc_illumina/desRot/46395_CTTGTA_L002_R2_fastqc.html",
        "outputs/fastqc_pacbio/macCal/macCal_filtered_subreads_fastqc.html",
        "outputs/fastqc_pacbio/desRot/desRot_filtered_subreads_fastqc.html"
    shell: """
        ssh {REMOTE_HOST} "mkdir -p {REMOTE_PATH}/bat"
        scp {input} {REMOTE_HOST}:{REMOTE_PATH}/bat
    """

rule publish_preqc:
    input:
        "outputs/sga/macCal/46394_GCCAAT_L001.pdf",
        "outputs/sga/macCal/46394_GCCAAT_L002.pdf",
        "outputs/sga/desRot/46395_CTTGTA_L001.pdf",
        "outputs/sga/desRot/46395_CTTGTA_L002.pdf",
        "outputs/sga/macCal/filtered_subreads.pdf",
        "outputs/sga/desRot/filtered_subreads.pdf"
    run:
        for i in input:
            path, fname = os.path.split(i)
            _, parent = os.path.split(path)
            shell('ssh {REMOTE_HOST} "mkdir -p {REMOTE_PATH}/bat/{parent}"')
            shell('rsync -e ssh {i} {REMOTE_HOST}:{REMOTE_PATH}/bat/{parent}/{fname}')

rule lordec_correct_all:
    input:
        "outputs/pacbio_lordec/macCal.fastq.gz",
        "outputs/pacbio_lordec/desRot.fastq.gz"

ruleorder: sga_index_pacbio > sga_index_illumina
