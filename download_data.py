#!/usr/bin/env python3
"""
Download UniProt Swiss-Prot human flat file for offline PTM queries.
The server serves the file gzip-compressed; this module handles decompression.
"""
import gzip
import os
import shutil
import requests
from datetime import datetime, timedelta

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
GZ_FILE = os.path.join(DATA_DIR, "uniprot_sprot_human.dat.gz")
FLAT_FILE = os.path.join(DATA_DIR, "uniprot_sprot_human.dat")
FLAT_FILE_URL = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
                 "knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat")
TTL_DAYS = 7
DOWNLOAD_TIMEOUT = 90
MIN_FLAT_FILE_SIZE = 100 * 1024 * 1024  # incomplete downloads must re-download


def needs_refresh():
    if not os.path.exists(FLAT_FILE):
        return True
    if os.path.getsize(FLAT_FILE) < MIN_FLAT_FILE_SIZE:
        return True  # interrupted/incomplete download
    mtime = datetime.fromtimestamp(os.path.getmtime(FLAT_FILE))
    return datetime.now() - mtime > timedelta(days=TTL_DAYS)


def _is_gzip(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def download(timeout=DOWNLOAD_TIMEOUT):
    os.makedirs(DATA_DIR, exist_ok=True)
    print("[*] Downloading UniProt human flat file (gzip, ~120 MB, may take minutes)...")
    r = requests.get(FLAT_FILE_URL, timeout=timeout, stream=True)
    r.raise_for_status()
    tmp = GZ_FILE + ".part"
    with open(tmp, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)

    if _is_gzip(tmp):
        print("[*] Decompressing...")
        with gzip.open(tmp, "rb") as gz, open(FLAT_FILE, "wb") as out:
            shutil.copyfileobj(gz, out)
        os.remove(tmp)
    else:
        os.replace(tmp, FLAT_FILE)

    print(f"[*] Offline data ready: {FLAT_FILE}")


def ensure_fresh(timeout=DOWNLOAD_TIMEOUT):
    """Best-effort refresh of offline data. Never raises; returns bool."""
    if not needs_refresh():
        return os.path.exists(FLAT_FILE)
    try:
        download(timeout=timeout)
        return True
    except Exception as e:
        print(f"[!] Offline data download failed: {e}")
        return False


def report_status():
    """Report offline data state without downloading. Returns bool (usable)."""
    if not os.path.exists(FLAT_FILE):
        print("[!] Offline UniProt flat file not found (use --refresh to download).")
        return False
    if os.path.getsize(FLAT_FILE) < MIN_FLAT_FILE_SIZE:
        print("[!] Offline UniProt flat file is incomplete (use --refresh to re-download).")
        return False
    return True
