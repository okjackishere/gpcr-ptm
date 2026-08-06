#!/usr/bin/env python3
"""
Update offline PTM data from dbPTM and iPTMnet.

Usage:
  python update_ptm_db.py --dbptm          # download dbPTM flat file
  python update_ptm_db.py --iptmnet        # download iPTMnet bulk ptm.txt
  python update_ptm_db.py --all            # download both

dbPTM URL is configurable via env var DBPTM_URL (defaults to the CUHK host).
"""
import argparse
import os

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
IPTMNET_PTM_FILE = os.path.join(DATA_DIR, "iptmnet_ptm.txt")
IPTMNET_PTM_URL = "https://research.bioinformatics.udel.edu/iptmnet_data/files/current/ptm.txt"


def download(url, path):
    import requests
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"[*] Downloading {url} ...")
    r = requests.get(url, timeout=180, stream=True)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
    print(f"[*] Saved to {path}")


def main():
    parser = argparse.ArgumentParser(description="Update offline PTM data (dbPTM / iPTMnet)")
    parser.add_argument("--dbptm", action="store_true", help="Download dbPTM flat file")
    parser.add_argument("--iptmnet", action="store_true", help="Download iPTMnet bulk ptm.txt")
    parser.add_argument("--all", action="store_true", help="Download both")
    parser.add_argument("--dbptm-url", default=os.environ.get("DBPTM_URL"), help="Override dbPTM URL")
    args = parser.parse_args()

    if args.all:
        args.dbptm = args.iptmnet = True

    if args.dbptm:
        from ptm_sources import download_dbptm
        try:
            download_dbptm(url=args.dbptm_url)
        except Exception as e:
            print(f"[!] dbPTM download failed: {e}")
            print("    Check the URL; override with --dbptm-url or env DBPTM_URL.")

    if args.iptmnet:
        try:
            download(IPTMNET_PTM_URL, IPTMNET_PTM_FILE)
        except Exception as e:
            print(f"[!] iPTMnet bulk download failed: {e}")


if __name__ == "__main__":
    main()
