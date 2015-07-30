WORKDIR = "/bat_data/leafnosed"

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

REMOTE_HOST = "athyra"
REMOTE_PATH = "public_html"
