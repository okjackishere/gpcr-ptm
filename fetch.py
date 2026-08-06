import requests

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


def fetch_uniprot(accession):
    url = f"{UNIPROT_BASE}/{accession}.json"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()


def fetch_iptmnet(accession):
    """Fetch iPTMnet PTM sites via the live REST API."""
    from ptm_sources import fetch_iptmnet_sites
    try:
        return fetch_iptmnet_sites(accession)
    except Exception:
        return []
