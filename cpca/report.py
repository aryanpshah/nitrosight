"""Render a CPCAResult as a human-readable, factor-by-factor report."""

from __future__ import annotations

import textwrap

from . import policy, rules
from .scorer import CPCAResult

RULE = "=" * 72
THIN = "-" * 72


def format_report(result: CPCAResult, name: str | None = None) -> str:
    lines: list[str] = [RULE]
    lines.append(f"  NDSRI POTENCY ASSESSMENT{'  -  ' + name if name else ''}")
    lines.append(RULE)
    lines.append(f"SMILES: {result.smiles}")
    lines.append("")

    if not result.applicable:
        lines.append("RESULT: CPCA does not apply")
        lines.append(f"Reason: {result.not_applicable_reason}")
        _append_warnings(lines, result)
        _append_footer(lines, result)
        return "\n".join(lines)

    lines.append(f"RESULT: Potency Category {result.category}"
                 f"    Recommended AI = {result.ai_ng_day:g} ng/day")
    lines.append(f"        {rules.CATEGORY_RATIONALE[result.category]}")
    lines.append("")

    multi = len(result.sites) > 1
    governing = result.governing_site

    for site in result.sites:
        header = f"N-nitroso site {site.site_id} (amine N at atom {site.amine_n_idx})"
        if multi and site is governing:
            header += "   <- governs the limit"
        lines.append(THIN)
        lines.append(header)
        lines.append(THIN)

        if site.alpha_h_counts:
            a, b = site.alpha_h_counts
            lines.append(f"  alpha-hydrogen counts (lowest first): {a},{b}")

        if site.triage_exit:
            lines.append(f"  Triage: {site.triage_exit}")
            lines.append(f"  -> Potency Category {site.category}, "
                         f"AI = {site.ai_ng_day:g} ng/day (no potency score needed)")
            lines.append("")
            continue

        if site.potency_score is None:
            for note in site.notes:
                lines.append(f"  ** {note}")
            lines.append("")
            continue

        lines.append("")
        lines.append("  Factors:")
        for f in site.factors:
            flag = "  [needs review]" if f.uncertain else ""
            label = f"{f.name}  [{f.reference}]{flag}"
            wrapped = textwrap.wrap(label, width=60) or [label]
            lines.append(f"    {f.score:+d}   {wrapped[0]}")
            for cont in wrapped[1:]:
                lines.append(f"         {cont}")
            if f.detail:
                for cont in textwrap.wrap(f.detail, width=60):
                    lines.append(f"           {cont}")

        lines.append("")
        calc = " + ".join(f"({f.score:+d})" for f in site.factors)
        lines.append(f"  Potency Score = {calc} = {site.potency_score}")
        lines.append(f"  -> Potency Category {site.category}, "
                     f"AI = {site.ai_ng_day:g} ng/day")
        lines.append("")

    _append_warnings(lines, result)

    if any(f.uncertain for s in result.sites for f in s.factors):
        lines.append(THIN)
        lines.append("ITEMS MARKED [needs review]")
        lines.append("  Electron-withdrawing group matching is the least precisely")
        lines.append("  specified part of CPCA. FDA defers to Cross & Ponting (2021)")
        lines.append("  rather than listing the groups. A qualified assessor should")
        lines.append("  confirm these matches before the result is relied on.")
        lines.append("")

    _append_footer(lines, result)
    return "\n".join(lines)


def _append_warnings(lines: list[str], result: CPCAResult) -> None:
    if result.warnings:
        lines.append(THIN)
        lines.append("WARNINGS")
        for w in result.warnings:
            lines.append(f"  ! {w}")
        lines.append("")


def _append_footer(lines: list[str], result: CPCAResult | None = None) -> None:
    lines.append(RULE)
    lines.append(f"Source: {rules.GUIDANCE_CITATION}")
    if result is not None:
        lines.append(_ruleset_line(result))
    lines.append("Screening aid for a qualified assessor. Not a regulatory determination.")
    lines.append(RULE)


def _ruleset_line(result: CPCAResult) -> str:
    """FDA's guidance and FDA's own reference code disagree on four rules, so every
    result has to say which reading produced it."""
    try:
        rs = policy.get(result.ruleset)
    except ValueError:
        return f"Ruleset: {result.ruleset}"

    followed = ", ".join(f"{rule} -> {source}" for rule, source in rs.divergences())
    return (f"Ruleset: {rs.name} ({rs.description})\n"
            f"  On the four guidance/reference-code divergences: {followed}")
