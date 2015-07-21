WORKDIR = "/mnt/research/ged/irberlui/biodata/bat"

with open(WORKDIR + "/files", 'r') as f:
    RAW_DATA = f.read().splitlines()

REMOTE_HOST = "athyra"
REMOTE_PATH = "public_html"
