import requests
import re

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def resolve_input(text):
    text = text.strip()

    # Looks like UniProt ID?
    if re.match(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$', text):
        return {"accession": text, "name": text, "gene": text}

    # Search by name/gene
    params = {
        "query": text,
        "fields": "accession,id,gene_names,protein_name,organism_id",
        "format": "json",
        "size": 50,
    }
    try:
        r = requests.get(UNIPROT_SEARCH, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"[!] UniProt search failed: {e}")
        return None

    results = data.get("results", [])
    if not results:
        return None

    # Filter for GPCRs by protein name
    gpcrs = []
    for res in results:
        pname = res.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "").lower()
        if "g-protein coupled" in pname or "receptor" in pname:
            gpcrs.append(res)

    candidates = gpcrs if gpcrs else results

    if len(candidates) == 1:
        c = candidates[0]
        genes = c.get("genes", [])
        gene = genes[0].get("geneName", {}).get("value", c.get("uniProtkbId", c["primaryAccession"])) if genes else c.get("uniProtkbId", c["primaryAccession"])
        return {
            "accession": c["primaryAccession"],
            "name": c.get("uniProtkbId", c["primaryAccession"]),
            "gene": gene,
        }

    # Multiple candidates: prefer human + reviewed
    def sort_key(c):
        is_reviewed = 1 if "reviewed" in c.get("entryType", "").lower() else 0
        is_human = 1 if c.get("organism", {}).get("taxonId") == 9606 else 0
        return (is_human, is_reviewed)

    candidates.sort(key=sort_key, reverse=True)
    best = candidates[0]
    genes = best.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", best.get("uniProtkbId", best["primaryAccession"])) if genes else best.get("uniProtkbId", best["primaryAccession"])
    pname = best.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    print(f"[*] Auto-selected: {best['primaryAccession']} ({gene}) — {pname}")
    return {
        "accession": best["primaryAccession"],
        "name": best.get("uniProtkbId", best["primaryAccession"]),
        "gene": gene,
    }
