"""
Which source of truth to follow where FDA's guidance and FDA's code disagree.

FDA published two things that do not match each other:

  * the guidance PDF, which describes the method in words
  * featurize-nitrosamines, the Java program that implements it

For four rules the two say different things, and there is no way to compute which is
right. Both are official. So the choice is made explicitly here, in one place, instead
of being buried in the scorer.

Each divergence is documented in docs/reference-implementation-divergences.md with the
guidance quote, the line of FDA source, and a worked compound.

Three presets:

  balanced (default)  Follow the PDF where it is clear and the code contradicts it.
                      Follow the code where the PDF is ambiguous and the code settles
                      it. This is the defensible reading: a reviewer checking us
                      against the published document finds us obeying it, and where
                      the document cannot be obeyed unambiguously we defer to the
                      operative tool.

  guidance            Follow the PDF everywhere it says anything. Use this to show
                      what the published document, read literally, produces.

  reference           Follow FDA's Java program everywhere. Use this to reproduce
                      featurize-nitrosamines, e.g. when benchmarking against it.

Not switchable: the electron-withdrawing group rule (Table B footnote **). The guidance
does not list the groups at all, it cites Cross & Ponting (2021), so there is no
"guidance" behaviour to follow. cpca/ewg.py reproduces the reference implementation's
ten patterns under every preset.
"""

from __future__ import annotations

from dataclasses import dataclass

GUIDANCE = "guidance"
REFERENCE = "reference"


@dataclass(frozen=True)
class Ruleset:
    """Which reading to apply at each of the four divergence points."""

    name: str
    description: str

    # Divergence 1: Table A has six entries in the PDF, ten in the code. The extra
    # 3,3 row is the one that matters; it makes the code score compounds the PDF puts
    # out of scope. GUIDANCE declines to score them, REFERENCE scores them.
    alpha_h_table: str = GUIDANCE

    # Divergence 2: Table A footnote $, the ethyl exception. The PDF requires it. The
    # code computes the flag and never uses it. True applies it.
    ethyl_exception: bool = True

    # Divergence 3: "tertiary alpha-carbon". The PDF says sp3 carbon bonded to three
    # other carbons. The code tests only for a carbon with no hydrogens and four
    # bonds, so a hydrogen-free alpha-carbon bearing a heteroatom qualifies.
    tertiary_alpha_carbon: str = REFERENCE

    # Divergence 4: beta-hydroxyl. The PDF requires an sp3 beta-carbon. The code
    # requires only that the beta-carbon carry no double bond, and an aromatic bond
    # is not a double bond, so a phenol counts under the code and not under the PDF.
    beta_hydroxyl: str = REFERENCE

    def divergences(self) -> list[tuple[str, str]]:
        """(rule, which source we followed), for the report and for audit trails."""
        return [
            ("Table A alpha-hydrogen patterns", self.alpha_h_table),
            ("Table A ethyl exception (footnote $)",
             GUIDANCE if self.ethyl_exception else REFERENCE),
            ("Tertiary alpha-carbon definition", self.tertiary_alpha_carbon),
            ("Beta-hydroxyl definition", self.beta_hydroxyl),
        ]


BALANCED = Ruleset(
    name="balanced",
    description=("Follow the guidance where it is clear and the code contradicts it; "
                 "follow the code where the guidance is ambiguous."),
    alpha_h_table=GUIDANCE,
    ethyl_exception=True,
    tertiary_alpha_carbon=REFERENCE,
    beta_hydroxyl=REFERENCE,
)

GUIDANCE_ONLY = Ruleset(
    name="guidance",
    description="Follow the published FDA guidance document throughout.",
    alpha_h_table=GUIDANCE,
    ethyl_exception=True,
    tertiary_alpha_carbon=GUIDANCE,
    beta_hydroxyl=GUIDANCE,
)

REFERENCE_ONLY = Ruleset(
    name="reference",
    description=("Reproduce FDA's featurize-nitrosamines reference implementation "
                 "throughout, including its departures from the guidance."),
    alpha_h_table=REFERENCE,
    ethyl_exception=False,
    tertiary_alpha_carbon=REFERENCE,
    beta_hydroxyl=REFERENCE,
)

DEFAULT = BALANCED
DEFAULT_RULESET = BALANCED  # exported under this name from the cpca package

PRESETS: dict[str, Ruleset] = {
    BALANCED.name: BALANCED,
    GUIDANCE_ONLY.name: GUIDANCE_ONLY,
    REFERENCE_ONLY.name: REFERENCE_ONLY,
}

PRESET_NAMES = list(PRESETS)


def get(ruleset: Ruleset | str | None) -> Ruleset:
    """Accept a Ruleset, a preset name, or None (meaning the default)."""
    if ruleset is None:
        return DEFAULT
    if isinstance(ruleset, Ruleset):
        return ruleset
    try:
        return PRESETS[ruleset]
    except KeyError:
        raise ValueError(
            f"Unknown ruleset {ruleset!r}. Choose one of: "
            f"{', '.join(PRESET_NAMES)}") from None
