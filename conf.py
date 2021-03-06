WORKDIR = "/mnt/research/ged/irberlui/biodata/bat"

with open(WORKDIR + "/files", 'r') as f:
    RAW_DATA = f.read().splitlines()

PACBIO_RAW = {
  'macCal': ['inputs/' + f for f in RAW_DATA if f.startswith("Pacbio/TK163824")],
  'desRot': ['inputs/' + f for f in RAW_DATA if f.startswith("Pacbio/TK169403")]
}

ILLUMINA_RAW = {
  'macCal': ['inputs/' + f for f in RAW_DATA if 'Sample_46394' in f],
  'desRot': ['inputs/' + f for f in RAW_DATA if 'Sample_46395' in f]
}

SAMPLE_TO_SPECIES = {
    '46394': 'macCal',
    '46395': 'desRot',
}

SPECIES_TO_SAMPLE = {v:k for (k, v) in SAMPLE_TO_SPECIES.items()}

SPECIES_TO_BARCODE = {
    'macCal': "GCCAAT",
    'desRot': "CTTGTA"
}



REMOTE_HOST = "athyra"
REMOTE_PATH = "public_html"
