#!/usr/bin/env python

"""Download a CSV from google sheets."""

import sys
import csv
from io import StringIO
from urllib.request import urlopen
import yaml

def download_google_sheet(name, path_yaml, path_out):
    """Download a CSV from google sheets."""
    with open(path_yaml, encoding="UTF8") as f_in:
        info = yaml.safe_load(f_in)
    for key in info:
        slug = info[key]["slug"]
        sheets = info[key]["sheets"]
        gid = sheets[name]
        url = "https://docs.google.com/spreadsheets/d/e/" \
            f"{slug}/pub?gid={gid}&single=true&output=csv"
        with urlopen(url) as f_in, open(path_out, "wt", encoding="UTF8") as f_out:
            reader = csv.reader(StringIO(f_in.read().decode("UTF8")))
            writer = csv.writer(f_out, lineterminator="\n")
            for row in reader:
                writer.writerow(row)

if __name__ == "__main__":
    download_google_sheet(sys.argv[1], sys.argv[2], sys.argv[3])
