"""
Confidence scoring for *predicted* (L3) PTM sites only.

score = base * (0.4 + 0.6 * conservation)

  - base:       motif specificity weight (see predict.py MOTIF_BASE)
  - conservation: fraction of orthologs sharing the residue at the aligned column

The score is only a relative ranking aid, NOT a probability. Evidence-backed
(L1/L2) sites are scored separately in verified.py and are never scored here.
The old "prior" term (same-type sites elsewhere in the protein boosting every
candidate) is removed: it was circular reasoning.

Confidence tier (explicit, not hidden in the score):
  High   base>=0.7 and conservation>=0.7
  Medium base>=0.55 and conservation>=0.5
  Low    otherwise
"""


def score_predictions(predictions, verified):
    for p in predictions:
        base = p.get("base", 0.3)
        cons = p.get("conservation")
        if cons is None:
            cons = 0.5
            if p.get("reasons"):
                p["reasons"] = p["reasons"] + ["保守性未知(无直系同源比对)"]
        score = round(min(1.0, base * (0.4 + 0.6 * cons)), 3)

        if base >= 0.7 and cons >= 0.7:
            confidence = "High"
        elif base >= 0.55 and cons >= 0.5:
            confidence = "Medium"
        else:
            confidence = "Low"

        p["score"] = score
        p["confidence"] = confidence

    predictions.sort(key=lambda x: x["score"], reverse=True)
    return predictions
