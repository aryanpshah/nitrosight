"""
Constants transcribed from FDA's Carcinogenic Potency Categorization Approach (CPCA).

Source: "Recommended Acceptable Intake Limits for Nitrosamine Drug Substance-Related
Impurities (NDSRIs)", FDA Guidance for Industry, August 2023.
Local copy: guidelines/FDA_nitrosamine_AI_limits.pdf

This module holds ONLY the guidance's own numbers and definitions, kept separate from
the algorithm so a reviewer can check it line-by-line against the PDF. Every entry
carries the figure/table it came from.
"""

from rdkit import Chem

# Figure 1 (p.10): potency category to recommended acceptable intake

CATEGORY_AI_NG_DAY: dict[int, float] = {
    1: 26.5,    # equal to the class-specific limit for nitrosamines, based on NDEA
    2: 100.0,   # representative of NDMA (96 ng/day) and NNK (100 ng/day)
    3: 400.0,   # 4-fold decrease in potency from Category 2
    4: 1500.0,  # TTC per ICH M7(R2)
    5: 1500.0,  # TTC per ICH M7(R2)
}

CATEGORY_RATIONALE: dict[int, str] = {
    1: "Predicted potency no higher than the class-specific nitrosamine limit (NDEA-based).",
    2: "Predicted potency no higher than NDMA and NNK.",
    3: "Lower potency than Category 2, e.g. a weakly deactivating feature is present.",
    4: "May be activated via alpha-hydroxylation but predicted to be of low potency.",
    5: "Not predicted to be metabolically activated via an alpha-hydroxylation pathway.",
}


def category_for_score(potency_score: int) -> int:
    """Figure 1, lower half. Note the guidance tests 'score <= 1', not 'score == 1',
    so activating features that drive the score below 1 still land in Category 1."""
    if potency_score >= 4:
        return 4
    if potency_score == 3:
        return 3
    if potency_score == 2:
        return 2
    return 1


# Table A (p.18-19): alpha-hydrogen score
# Key is the count of hydrogens on each alpha-carbon, LOWEST FIRST.

ALPHA_HYDROGEN_SCORE: dict[tuple[int, int], int] = {
    (0, 2): 3,  # see ETHYL_EXCEPTION below
    (0, 3): 2,
    (1, 2): 3,
    (1, 3): 3,
    (2, 2): 1,
    (2, 3): 1,
}

# Table A footnote $: "A score of 3 applies when the methylene alpha-carbon is not part
# of an ethyl group. If the methylene alpha-carbon is part of an ethyl group, a score of
# 2 should be applied."  Applies to the (0,2) row.
ETHYL_EXCEPTION_PAIR = (0, 2)
ETHYL_EXCEPTION_SCORE = 2

# Combinations deliberately absent from Table A, (3,3) among them. The guidance says a
# structure "is expected to match exactly one of the alpha-hydrogen definitions in Table
# A", so anything else is out of scope and must be reported, never guessed.

# FDA's reference implementation (featurize-nitrosamines,
# FeaturizeNitrosamine.java line 733) builds a TEN-entry lookup instead of Table A's
# six. The four extra rows are not in the published guidance:
#
#   (0,0) -> 5, (0,1) -> 4, (1,1) -> 4     unreachable here: Figure 1's triage gates
#                                          catch all three before Table A is consulted
#   (3,3) -> 1                             reachable, and the one that matters: it
#                                          makes the code score N-nitrosodimethylamine
#                                          and its relatives, which the guidance text
#                                          places outside the method's scope
#
# Selected by Ruleset.alpha_h_table == "reference". See cpca/policy.py and
# docs/reference-implementation-divergences.md.
REFERENCE_ALPHA_HYDROGEN_SCORE: dict[tuple[int, int], int] = {
    **ALPHA_HYDROGEN_SCORE,
    (0, 0): 5,
    (0, 1): 4,
    (1, 1): 4,
    (3, 3): 1,
}


# Table B (p.19-20): deactivating features (positive scores)

CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX1-]")

RING_FEATURE_PYRROLIDINE = 3      # N-nitroso group in a pyrrolidine ring
RING_FEATURE_6RING_SULFUR = 3     # 6-membered ring containing at least one sulfur
RING_FEATURE_5_OR_6_RING = 2      # any other 5- or 6-membered ring
RING_FEATURE_MORPHOLINE = 1       # morpholine ring
RING_FEATURE_7RING = 1            # 7-membered ring

CARBOXYLIC_ACID_SCORE = 3
LONG_CHAIN_BOTH_SIDES_SCORE = 1   # chains of >=5 heavy atoms on both sides, acyclic N-nitroso
CHAIN_MIN_LENGTH = 5
CHAIN_MAX_ATOMS_IN_ONE_RING = 4

EWG_ONE_SIDE_SCORE = 1
EWG_BOTH_SIDES_SCORE = 2

BETA_HYDROXYL_ONE_SIDE_SCORE = 1
BETA_HYDROXYL_BOTH_SIDES_SCORE = 2

# Table B footnote **: electron-withdrawing groups EXCLUDE carboxylic acid and aryl
# (counted separately) and ketone (conflicting data). The guidance defers to
# Cross KP & Ponting DJ (2021), Comput Toxicol 20:100186, for the full list.
#
# CAUTION: this is the single least precisely specified part of CPCA. The list below is
# a deliberate, conservative reading; the scorer flags every structure where it fires so
# a qualified assessor can confirm. Do not treat these matches as settled.
ELECTRON_WITHDRAWING_GROUPS: list[tuple[str, str]] = [
    ("ester (C-bonded)",   "[CX3](=[OX1])[OX2][#6]"),
    ("ester (O-bonded)",   "[OX2][CX3]=[OX1]"),
    ("amide",              "[CX3](=[OX1])[NX3]"),
    ("nitrile",            "[NX1]#[CX2]"),
    ("nitro",              "[NX3+](=[OX1])[OX1-]"),
    ("sulfonyl / sulfone", "[SX4](=[OX1])(=[OX1])"),
    ("sulfoxide",          "[SX3](=[OX1])"),
    ("halogen",            "[F,Cl,Br,I]"),
    ("trifluoromethyl",    "[CX4](F)(F)F"),
]

EWG_PATTERNS = [(name, Chem.MolFromSmarts(sma)) for name, sma in ELECTRON_WITHDRAWING_GROUPS]


# Table C (p.20): activating features (negative scores)

ARYL_ON_ALPHA_CARBON_SCORE = -1   # benzylic / pseudo-benzylic substituent
METHYL_ON_BETA_CARBON_SCORE = -1


# Applicability (section IV.C, p.9)

# "applies to NDSRIs bearing a carbon atom on both sides of the N-nitroso group, where
#  the carbon is not directly double bonded to a heteroatom (i.e., N-nitrosamides,
#  N-nitrosoureas, N-nitrosoguanidines and other related structures are excluded).
#  Additionally, [it] does not apply to NDSRIs where the N-nitroso group is within an
#  aromatic ring."
NITROSAMINE_CORE = Chem.MolFromSmarts("[NX3](-[#6])(-[#6])-[NX2]=[OX1]")

GUIDANCE_CITATION = (
    "FDA, Recommended Acceptable Intake Limits for Nitrosamine Drug Substance-Related "
    "Impurities (NDSRIs), Guidance for Industry, August 2023"
)
