"""
End-to-end assessment of one impurity.

Shared by the CLI and the desktop app so both give identical answers:

    input (SMILES / name / CAS)
        -> structure
        -> CPCA prediction
        -> FDA published limit, if one exists, which overrides the prediction
        -> optional ppm control limit for a given maximum daily dose
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from rdkit import Chem, RDLogger

from cpca import CPCAResult, format_report, score_smiles
from resolve import FDALimits, PubChemClient, PublishedLimit

RDLogger.DisableLog("rdApp.*")

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


@dataclass
class ScanResult:
    query: str
    input_kind: str = "unknown"          # smiles | name | cas
    smiles: str | None = None
    resolved_name: str | None = None
    pubchem_cid: int | None = None
    cpca: CPCAResult | None = None
    published: PublishedLimit | None = None
    error: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def final_ai_ng_day(self) -> float | None:
        if self.published is not None:
            return self.published.ai_ng_day
        if self.cpca is not None and self.cpca.applicable:
            return self.cpca.ai_ng_day
        return None

    @property
    def final_category(self) -> int | None:
        if self.published is not None and self.published.category is not None:
            return self.published.category
        if self.cpca is not None and self.cpca.applicable:
            return self.cpca.category
        return None

    @property
    def final_basis(self) -> str:
        if self.published is not None:
            return {"published": "FDA published limit",
                    "surrogate": "FDA limit (read-across from a surrogate)",
                    "interim": "FDA interim limit (time-limited)"}.get(
                        self.published.basis, "FDA published limit")
        if self.cpca is not None and self.cpca.applicable:
            return "CPCA prediction"
        return "no limit determined"

    @property
    def overridden(self) -> bool:
        """Did a published limit displace a different CPCA prediction?"""
        return (self.published is not None
                and self.cpca is not None and self.cpca.applicable
                and self.published.category is not None
                and self.published.category != self.cpca.category)

    def control_limit_ppm(self, max_daily_dose_mg: float) -> float | None:
        """ppm = AI (ng/day) / MDD (mg). FDA, Control of Nitrosamine Impurities."""
        ai = self.final_ai_ng_day
        if ai is None or max_daily_dose_mg <= 0:
            return None
        return ai / max_daily_dose_mg


def classify_input(text: str) -> str:
    text = text.strip()
    if CAS_RE.match(text):
        return "cas"
    if Chem.MolFromSmiles(text) is not None:
        return "smiles"
    return "name"


def scan(query: str, client: PubChemClient | None = None,
         limits: FDALimits | None = None, kind: str = "auto",
         ruleset=None) -> ScanResult:
    """Resolve, score, and apply any published limit. Network access only when the
    input is a name or CAS number.

    `ruleset` is passed straight to the scorer: it selects what to follow where FDA's
    guidance and FDA's reference code disagree. See cpca/policy.py."""
    query = query.strip()
    if not query:
        return ScanResult(query=query, error="Nothing to score.")

    kind = classify_input(query) if kind == "auto" else kind
    result = ScanResult(query=query, input_kind=kind)

    if kind == "smiles":
        result.smiles = query
    else:
        client = client or PubChemClient(verbose=False)
        hit = client.lookup(query)
        client.save()
        if not hit.found:
            result.error = (
                f"PubChem has no match for {query!r}. "
                "Try a CAS number, the parent drug name, or paste a SMILES string.")
            return result
        result.smiles = hit.smiles
        result.resolved_name = hit.iupac_name
        result.pubchem_cid = hit.cid

    if Chem.MolFromSmiles(result.smiles) is None:
        result.error = "That structure could not be parsed."
        return result

    result.cpca = score_smiles(result.smiles, ruleset)

    # A published limit governs; the prediction is kept for comparison.
    limits = limits or FDALimits()
    if stale := limits.staleness_warning():
        result.notes.append(stale)
    for candidate in filter(None, (query if kind == "name" else None,
                                   result.resolved_name)):
        if (pub := limits.lookup(candidate)) is not None:
            result.published = pub
            break

    return result


def format_scan_report(result: ScanResult | None,
                       max_daily_dose_mg: float | None = None) -> str:
    """A saveable text report: the governing limit, then the full CPCA breakdown."""
    if result is None:
        return "Nothing assessed."

    rule = "=" * 72
    lines = [rule, "  NITROSAMINE IMPURITY ASSESSMENT", rule,
             f"Query        : {result.query}",
             f"Read as      : {result.input_kind}"]
    if result.smiles:
        lines.append(f"SMILES       : {result.smiles}")
    if result.resolved_name:
        lines.append(f"Resolved name: {result.resolved_name}")
    if result.pubchem_cid:
        lines.append(f"PubChem CID  : {result.pubchem_cid}")
    lines.append("")

    if result.error:
        lines += [f"NOT ASSESSED: {result.error}", rule]
        return "\n".join(lines)

    ai = result.final_ai_ng_day
    if ai is None:
        reason = (result.cpca.not_applicable_reason if result.cpca
                  else "no limit could be determined")
        lines += [f"NO LIMIT DETERMINED: {reason}", rule]
        return "\n".join(lines)

    cat = result.final_category
    lines.append(f"GOVERNING LIMIT : {ai:g} ng/day"
                 + (f"   (Potency Category {cat})" if cat else ""))
    lines.append(f"BASIS           : {result.final_basis}")

    if result.published is not None:
        lines.append(f"FDA entry       : {result.published.name}")
        if result.published.source:
            lines.append(f"  source        : {result.published.source}")
        if result.published.surrogate:
            lines.append(f"  surrogate     : {result.published.surrogate}")
        if result.published.is_interim:
            lines.append("  NOTE          : INTERIM limit"
                         + (f" ({result.published.duration})"
                            if result.published.duration else ""))
    if result.overridden and result.cpca:
        lines.append(f"  supersedes    : CPCA prediction of Category "
                     f"{result.cpca.category} ({result.cpca.ai_ng_day:g} ng/day)")

    if max_daily_dose_mg:
        ppm = result.control_limit_ppm(max_daily_dose_mg)
        lines += ["",
                  f"CONTROL LIMIT   : {ppm:.4g} ppm  "
                  f"at a max daily dose of {max_daily_dose_mg:g} mg",
                  "  ppm = AI (ng/day) / MDD (mg)"]

    for note in result.notes:
        lines.append(f"! {note}")

    lines += ["", rule, "  CPCA BREAKDOWN", ""]
    lines.append(format_report(result.cpca, result.query) if result.cpca
                 else "no CPCA result")
    return "\n".join(lines)
