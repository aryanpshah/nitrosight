"""
The ruleset switch: does each of the four guidance/reference-code divergences
actually change behaviour, and does the default behave as documented?

Every worked compound here comes from docs/reference-implementation-divergences.md,
so this file is also the executable form of that document.

    python test_policy.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cpca import PRESET_NAMES, score_smiles                      # noqa: E402
from cpca.policy import BALANCED, DEFAULT, GUIDANCE, REFERENCE   # noqa: E402
from cpca.policy import GUIDANCE_ONLY, REFERENCE_ONLY, Ruleset, get  # noqa: E402

PASS, FAIL = 0, 0


def check(label: str, got, want) -> None:
    global PASS, FAIL
    if got == want:
        PASS += 1
        print(f"  [PASS] {label}")
    else:
        FAIL += 1
        print(f"  [FAIL] {label}\n         got {got!r}, want {want!r}")


def category(smiles: str, ruleset) -> int | None:
    """Category, or None when the scorer declines to score the structure."""
    r = score_smiles(smiles, ruleset)
    return r.category if r.applicable else None


print("=" * 74)
print("  THE DEFAULT IS THE DOCUMENTED ONE")
print("=" * 74)

check("default ruleset is 'balanced'", DEFAULT.name, "balanced")
check("no ruleset given means the default",
      get(None).name, DEFAULT.name)
check("a preset name resolves", get("reference").name, "reference")
check("results record which ruleset produced them",
      score_smiles("CCN(N=O)CC").ruleset, "balanced")

# The rule the default encodes: guidance where it is clear and the code contradicts
# it, reference code where the guidance is ambiguous.
check("balanced follows the guidance on Table A patterns",
      BALANCED.alpha_h_table, GUIDANCE)
check("balanced applies the guidance's ethyl exception",
      BALANCED.ethyl_exception, True)
check("balanced follows the reference code on tertiary alpha-carbon",
      BALANCED.tertiary_alpha_carbon, REFERENCE)
check("balanced follows the reference code on beta-hydroxyl",
      BALANCED.beta_hydroxyl, REFERENCE)

try:
    get("nonsense")
    check("an unknown ruleset name is rejected", "no error", "ValueError")
except ValueError as exc:
    check("an unknown ruleset name is rejected",
          "reference" in str(exc) and "nonsense" in str(exc), True)


print()
print("=" * 74)
print("  DIVERGENCE 1: Table A has six entries in the PDF, ten in the code")
print("=" * 74)

# N-nitrosodimethylamine. Alpha-hydrogen pattern 3,3, which Table A does not list.
NDMA = "CN(C)N=O"
check("guidance declines to score the 3,3 pattern",
      category(NDMA, GUIDANCE_ONLY), None)
check("reference code scores 3,3 as Category 1",
      category(NDMA, REFERENCE_ONLY), 1)
check("the default declines, following the guidance",
      category(NDMA, None), None)


print()
print("=" * 74)
print("  DIVERGENCE 2: the ethyl exception, Table A footnote $")
print("=" * 74)

# N-ethyl-N-(o-tolyl)nitrous amide. Pattern 0,2 with the methylene in an ethyl group,
# so footnote $ lowers the alpha-hydrogen score from 3 to 2.
ETHYL = "CCN(C1=CC=CC=C1C)N=O"
check("guidance applies footnote $ -> Category 2 (100 ng/day)",
      category(ETHYL, GUIDANCE_ONLY), 2)
check("reference code omits footnote $ -> Category 3 (400 ng/day)",
      category(ETHYL, REFERENCE_ONLY), 3)
check("the default applies it, following the guidance",
      category(ETHYL, None), 2)
check("the four-fold difference in the limit is real",
      (score_smiles(ETHYL, "guidance").ai_ng_day,
       score_smiles(ETHYL, "reference").ai_ng_day), (100.0, 400.0))


print()
print("=" * 74)
print("  DIVERGENCE 3: what counts as a tertiary alpha-carbon")
print("=" * 74)

# N-nitroso-2,2,5-trimethyl-1,3-oxazolidine. The alpha-carbon has no hydrogens but
# only two carbon substituents, the third being the ring oxygen.
OXAZ = "CC1CN(C(O1)(C)C)N=O"
check("guidance: not tertiary, so the compound is scored -> Category 4",
      category(OXAZ, GUIDANCE_ONLY), 4)
check("reference code: tertiary, so triage sends it to Category 5",
      category(OXAZ, REFERENCE_ONLY), 5)
check("the default follows the reference code",
      category(OXAZ, None), 5)

# A genuinely tertiary alpha-carbon reads the same either way.
TBUTYL = "CC(C)(C)N(C)N=O"
check("an all-carbon tertiary alpha-carbon is tertiary under both readings",
      (category(TBUTYL, GUIDANCE_ONLY), category(TBUTYL, REFERENCE_ONLY)), (5, 5))


print()
print("=" * 74)
print("  DIVERGENCE 4: beta-hydroxyl, where the code counts phenols")
print("=" * 74)

# A terminal OH on the beta-carbon satisfies both readings.
TERMINAL_OH = "OCCN(N=O)CCO"
guid = score_smiles(TERMINAL_OH, "guidance")
ref = score_smiles(TERMINAL_OH, "reference")


def has_beta_oh(result) -> bool:
    site = result.governing_site
    return bool(site) and any("Hydroxyl group bonded to beta-carbon" in f.name
                              for f in site.factors)


check("a terminal beta-OH is found under the guidance reading",
      has_beta_oh(guid), True)
check("a terminal beta-OH is found under the reference reading",
      has_beta_oh(ref), True)
check("both readings agree on the category there",
      guid.category == ref.category, True)

# A phenol: the hydroxyl sits on an AROMATIC beta-carbon. The guidance requires the
# beta-carbon to be sp3, which excludes it. The reference code tests only that the
# beta-carbon carries no double bond, and an aromatic bond is not a double bond, so
# it counts. This is the whole of divergence 4, and phenols are common in drugs.
PHENOL = "CN(N=O)c1ccccc1O"
check("a phenolic beta-OH is counted only by the reference reading",
      (has_beta_oh(score_smiles(PHENOL, "guidance")),
       has_beta_oh(score_smiles(PHENOL, "reference"))),
      (False, True))
check("and it changes the category",
      (category(PHENOL, GUIDANCE_ONLY), category(PHENOL, REFERENCE_ONLY)), (2, 3))

# An sp3 beta-carbon with a secondary hydroxyl reads the same either way, so the
# divergence really is about aromaticity and not about where the OH sits.
SECONDARY_OH = "CC(O)CN(N=O)CC"
check("a secondary beta-OH on an sp3 carbon is counted under both readings",
      (has_beta_oh(score_smiles(SECONDARY_OH, "guidance")),
       has_beta_oh(score_smiles(SECONDARY_OH, "reference"))),
      (True, True))


print()
print("=" * 74)
print("  ISOTOPE LABELS DO NOT CHANGE THE SCORE")
print("=" * 74)

# RDKit holds [2H] as an explicit atom, so a CD3 group reports zero hydrogens and
# Figure 1's first gate would send these to Category 5, 1500 ng/day. Each pair below
# is a labelled standard and its unlabelled parent; they must agree.
PAIRS = [
    ("N-nitrosoethylmethylamine-D3",
     "[2H]C([2H])([2H])N(CC)N=O", "CCN(C)N=O"),
    ("N-nitrosodipropylamine-D14",
     "[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])N(C([2H])([2H])C([2H])([2H])"
     "C([2H])([2H])[2H])N=O", "CCCN(CCC)N=O"),
    ("N-nitrosomorpholine-D8",
     "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])N(N=O)C1([2H])[2H]", "O=NN1CCOCC1"),
    ("N-nitrososarcosine-D3",
     "[2H]C([2H])([2H])N(CC(=O)O)N=O", "CN(CC(=O)O)N=O"),
]

for label, labelled, plain in PAIRS:
    check(f"{label} scores as its unlabelled parent",
          category(labelled, None), category(plain, None))

check("N-nitrosoethylmethylamine-D3 is Category 1, not Category 5",
      category("[2H]C([2H])([2H])N(CC)N=O", None), 1)

iso = score_smiles("[2H]C([2H])([2H])N(CC)N=O")
check("the report says the labels were removed",
      any("Isotope labels removed" in w for w in iso.warnings), True)
check("an unlabelled structure gets no such warning",
      any("Isotope" in w for w in score_smiles("CCN(C)N=O").warnings), False)

# A carbon-13 label changes no hydrogen count, so it must change no score either.
check("a 13C label is harmless",
      category("[13CH3]N(C)N=O", "reference"), category("CN(C)N=O", "reference"))

# The fix must not disturb salt stripping, which runs just before it.
check("isotope stripping coexists with salt stripping",
      category("[2H]C([2H])([2H])N(CC)N=O.Cl", None), 1)


print()
print("=" * 74)
print("  THE COUNTER-ION IS NEVER SCORED IN PLACE OF THE NITROSAMINE")
print("=" * 74)

NDEA = "CCN(N=O)CC"
check("plain NDEA is Category 1", category(NDEA, None), 1)

# Succinic acid has MORE heavy atoms than NDEA, so picking the largest fragment would
# score the counter-ion and report that CPCA does not apply.
BIG_SALT = "CCN(N=O)CC.OC(=O)CCC(=O)O"
check("a counter-ion larger than the nitrosamine is still stripped",
      category(BIG_SALT, None), 1)
check("and the carboxylic acid on it does not reach Table B",
      score_smiles(BIG_SALT).governing_site.potency_score,
      score_smiles(NDEA).governing_site.potency_score)
check("the report explains the fragment choice",
      any("NOT the nitrosamine" in w for w in score_smiles(BIG_SALT).warnings), True)

# A small counter-ion must behave exactly as before.
check("a small counter-ion is unaffected", category("CCN(N=O)CC.Cl", None), 1)

# When no fragment carries an N-nitroso group, fall back to the largest and decline.
check("a structure with no nitrosamine anywhere still declines",
      category("OC(=O)CCC(=O)O.Cl", None), None)


print()
print("=" * 74)
print("  AN N-NITROSO GROUP ON AN AROMATIC RING NITROGEN IS RECOGNISED")
print("=" * 74)

# The CPCA does not apply to these, but the tool has to find the group before it
# can say so. Detecting it needs BOTH [#7X3] instead of [NX3] and "~" instead of
# "-", because an aromatic ring nitrogen is not aliphatic and its ring bonds are
# not single bonds. Fixing either alone matches nothing, so both are tested.
AROMATIC_RING_N = {
    "N-nitroso esomeprazole": "COc1ccc2nc([S@](=O)Cc3ncc(C)c(OC)c3C)n(N=O)c2c1",
    "N-nitroso imidazole": "O=Nn1ccnc1",
    "N-nitroso benzimidazole": "O=NN1C=Nc2ccccc21",
}

for label, smi in AROMATIC_RING_N.items():
    res = score_smiles(smi)
    check(f"{label} is still out of scope", res.applicable, False)
    check(f"{label} is told the real reason",
          "aromatic ring" in res.not_applicable_reason.lower(), True)
    check(f"{label} is not told the group is missing",
          "no n-nitrosamine group found" in res.not_applicable_reason.lower(), False)

# An exocyclic nitrogen carrying two aromatic carbons is a different case and was
# always detected. It must keep scoring.
check("N-nitrosodiphenylamine still scores",
      category("O=NN(c1ccccc1)c1ccccc1", None), 5)
check("NDEA is unaffected by the broader pattern",
      category("CCN(N=O)CC", None), 1)


print()
print("=" * 74)
print("  THE SWITCH IS USABLE")
print("=" * 74)

check("three presets are exposed", sorted(PRESET_NAMES),
      ["balanced", "guidance", "reference"])

# A caller can mix and match rather than taking a preset wholesale.
custom = Ruleset(name="custom", description="guidance Table A, reference elsewhere",
                 alpha_h_table=GUIDANCE, ethyl_exception=False,
                 tertiary_alpha_carbon=REFERENCE, beta_hydroxyl=REFERENCE)
check("a hand-built ruleset works", category(ETHYL, custom), 3)
check("and is recorded by name", score_smiles(ETHYL, custom).ruleset, "custom")

check("every divergence is named in the audit trail",
      len(BALANCED.divergences()), 4)

# A compound none of the divergences touch must read the same under all three.
NDEA = "CCN(N=O)CC"
check("NDEA is unaffected by the choice",
      len({category(NDEA, p) for p in PRESET_NAMES}), 1)


print()
print("=" * 74)
print(f"  {PASS}/{PASS + FAIL} passed")
print("=" * 74)
raise SystemExit(1 if FAIL else 0)
