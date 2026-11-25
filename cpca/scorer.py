"""
FDA CPCA potency scorer.

Takes a SMILES string for a nitrosamine drug substance-related impurity (NDSRI) and
returns its predicted carcinogenic potency category, the recommended acceptable intake
limit, and every individual factor that contributed to the score.

The scorer is deterministic: the same structure always produces the same numbers, and
each number cites the table it came from. No language model is involved in the
calculation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem

from . import ewg, policy, rules
from .policy import GUIDANCE, Ruleset


# Result types


@dataclass
class Factor:
    """One line item in the potency score."""

    name: str
    score: int
    reference: str
    detail: str = ""
    uncertain: bool = False

    def __str__(self) -> str:
        flag = "  [needs review]" if self.uncertain else ""
        return f"{self.score:+d}  {self.name} ({self.reference}){flag}"


@dataclass
class SiteResult:
    """CPCA outcome for a single N-nitroso group."""

    site_id: int
    amine_n_idx: int
    alpha_h_counts: tuple[int, int] | None = None
    alpha_h_score: int | None = None
    # The three Figure 1 triage answers, recorded individually so they can be compared
    # against a manual assessment. All three are evaluated even though Figure 1 exits
    # early, so a hand-scored sheet can be checked question by question.
    gate_any_alpha_h: bool | None = None
    gate_multiple_alpha_h: bool | None = None
    gate_tertiary_alpha_c: bool | None = None
    factors: list[Factor] = field(default_factory=list)
    potency_score: int | None = None
    category: int = 0
    ai_ng_day: float = 0.0
    triage_exit: str | None = None
    notes: list[str] = field(default_factory=list)


@dataclass
class CPCAResult:
    """CPCA outcome for a whole molecule."""

    smiles: str
    applicable: bool
    category: int | None = None
    ai_ng_day: float | None = None
    sites: list[SiteResult] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    not_applicable_reason: str = ""
    ruleset: str = policy.DEFAULT.name

    @property
    def governing_site(self) -> SiteResult | None:
        """The N-nitroso group that determines the molecule's limit."""
        scored = [s for s in self.sites if s.category]
        return min(scored, key=lambda s: s.category) if scored else None


# Structure helpers


def _alpha_carbons(mol: Chem.Mol, amine_n_idx: int, nitroso_n_idx: int) -> list[Chem.Atom]:
    n = mol.GetAtomWithIdx(amine_n_idx)
    return [a for a in n.GetNeighbors()
            if a.GetAtomicNum() == 6 and a.GetIdx() != nitroso_n_idx]


def _beta_carbons(alpha: Chem.Atom, amine_n_idx: int) -> list[Chem.Atom]:
    return [a for a in alpha.GetNeighbors()
            if a.GetAtomicNum() == 6 and a.GetIdx() != amine_n_idx]


def _is_tertiary_alpha_carbon(atom: Chem.Atom, ruleset: Ruleset) -> bool:
    """Figure 1, gate 3. Divergence 3 of 4 - see cpca/policy.py.

    The guidance footnote defines this as "an alpha-carbon atom in an sp3 hybridization
    state, bonded to three other carbon atoms". FDA's reference implementation does NOT
    check that the three substituents are carbons - it tests only for a carbon with
    four bonds carrying no hydrogens:

        atn.isCarbon() && atn.getHCount() == 0 && atn.getBonds().count() == 4

    The two readings part company whenever a hydrogen-free alpha-carbon bears a
    heteroatom, as in N-nitroso-2,2,5-trimethyl-1,3-oxazolidine, where the alpha-carbon
    carries a ring oxygen and two methyls. The guidance text says that is not a tertiary
    alpha-carbon; the reference implementation says it is, and sends it to Category 5.
    """
    if atom.GetAtomicNum() != 6 or atom.GetTotalNumHs() != 0:
        return False

    if ruleset.tertiary_alpha_carbon == GUIDANCE:
        return (atom.GetHybridization() == Chem.HybridizationType.SP3
                and sum(1 for nb in atom.GetNeighbors()
                        if nb.GetAtomicNum() == 6) == 3)

    return atom.GetDegree() == 4


def _is_methyl(atom: Chem.Atom) -> bool:
    return (atom.GetAtomicNum() == 6
            and atom.GetTotalNumHs() == 3
            and atom.GetDegree() == 1)


def _in_ethyl_group(alpha: Chem.Atom, amine_n_idx: int) -> bool:
    """True when this methylene alpha-carbon is the CH2 of a terminal ethyl (-CH2-CH3)."""
    if alpha.GetTotalNumHs() != 2:
        return False
    others = _beta_carbons(alpha, amine_n_idx)
    return len(others) == 1 and _is_methyl(others[0])


def _has_double_bond(atom: Chem.Atom) -> bool:
    return any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())


def _has_hydroxyl(atom: Chem.Atom) -> bool:
    return any(nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() == 1 and nb.GetDegree() == 1
               for nb in atom.GetNeighbors())


def _is_terminal_heteroatom(atom: Chem.Atom) -> bool:
    """FDA AtomTest.isTerminalHeteroAtom: not carbon, at most one heavy neighbour."""
    if atom.GetAtomicNum() == 6:
        return False
    return sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1) <= 1


def _methyl_on_beta_carbon(amine_n_idx: int, alphas: list[Chem.Atom]) -> bool:
    """Table C: "Methyl group bonded to beta-carbon (cyclic or acyclic)", -1.

    The guidance wording is ambiguous about what counts. FDA's reference
    implementation (featurize-nitrosamines, FeaturizeNitrosamine.calculateFeatures)
    requires the beta-carbon to be a METHINE - exactly one hydrogen - so the methyl has
    to be a branch, not the far end of a straight chain. Reading it the looser way
    counts the terminal methyl of every ethyl and propyl group, which disagreed with
    FDA on 29 of our catalogue compounds.

    Their chain of filters, reproduced:
      alpha-carbon sp3, and carrying only carbon/hydrogen besides the nitrogen
      -> beta-carbon sp3, exactly one hydrogen, no terminal heteroatom attached
      -> does it bear a methyl?
    """
    for alpha in alphas:
        if alpha.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        if any(nb.GetIdx() != amine_n_idx and nb.GetAtomicNum() not in (1, 6)
               for nb in alpha.GetNeighbors()):
            continue

        for beta in alpha.GetNeighbors():
            if beta.GetAtomicNum() != 6 or beta.GetIdx() == amine_n_idx:
                continue
            if beta.GetHybridization() != Chem.HybridizationType.SP3:
                continue
            if beta.GetTotalNumHs() != 1:
                continue
            if any(_is_terminal_heteroatom(nb) for nb in beta.GetNeighbors()):
                continue
            if any(_is_methyl(nb) for nb in beta.GetNeighbors()
                   if nb.GetIdx() != alpha.GetIdx()):
                return True
    return False


def _chain_of_length_exists(mol: Chem.Mol, start: Chem.Atom, blocked_idx: int,
                            min_len: int, max_in_one_ring: int) -> bool:
    """Is there a simple path of `min_len` consecutive heavy atoms starting at `start`
    and walking away from `blocked_idx`, with no more than `max_in_one_ring` of those
    atoms sharing a single ring? (Table B, final row.)"""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]

    def ring_ok(path: list[int]) -> bool:
        pset = set(path)
        return all(len(pset & ring) <= max_in_one_ring for ring in rings)

    def walk(atom: Chem.Atom, path: list[int]) -> bool:
        if len(path) >= min_len:
            return ring_ok(path)
        for nb in atom.GetNeighbors():
            idx = nb.GetIdx()
            if idx == blocked_idx or idx in path or nb.GetAtomicNum() == 1:
                continue
            if walk(nb, path + [idx]):
                return True
        return False

    return walk(start, [start.GetIdx()])


def _classify_ring(mol: Chem.Mol, amine_n_idx: int) -> tuple[int, str] | None:
    """Score the smallest ring containing the N-nitroso nitrogen (Table B preamble:
    'the feature score for only the smallest matching ring should be applied')."""
    rings = [r for r in mol.GetRingInfo().AtomRings() if amine_n_idx in r]
    if not rings:
        return None

    ring = min(rings, key=len)
    others = [mol.GetAtomWithIdx(i) for i in ring if i != amine_n_idx]
    size = len(ring)

    if size == 5:
        if all(a.GetAtomicNum() == 6 and not a.GetIsAromatic() for a in others):
            return rules.RING_FEATURE_PYRROLIDINE, "N-nitroso group in a pyrrolidine ring"
        return rules.RING_FEATURE_5_OR_6_RING, "N-nitroso group in a 5-membered ring"

    if size == 6:
        if any(a.GetAtomicNum() == 16 for a in others):
            return (rules.RING_FEATURE_6RING_SULFUR,
                    "N-nitroso group in a 6-membered ring containing sulfur")
        oxygens = sum(1 for a in others if a.GetAtomicNum() == 8)
        carbons = sum(1 for a in others if a.GetAtomicNum() == 6)
        if oxygens == 1 and carbons == 4:
            return rules.RING_FEATURE_MORPHOLINE, "N-nitroso group in a morpholine ring"
        return rules.RING_FEATURE_5_OR_6_RING, "N-nitroso group in a 6-membered ring"

    if size == 7:
        return rules.RING_FEATURE_7RING, "N-nitroso group in a 7-membered ring"

    return None


# Feature scoring


def _beta_hydroxyl_sides(amine_n_idx: int, alphas: list[Chem.Atom],
                         ruleset: Ruleset) -> int:
    """How many sides of the N-nitroso group carry a beta-hydroxyl (Table B).
    Divergence 4 of 4 - see cpca/policy.py.

    The guidance requires the hydroxyl to sit on an sp3 beta-carbon. FDA's reference
    implementation instead asks only that the beta-carbon carry no double bond, and
    an aromatic bond is not a double bond. So a PHENOL - a hydroxyl on an aromatic
    beta-carbon - counts under the reference reading and not under the guidance's.
    That single case is the whole of this divergence, and it matters because phenol
    groups are common in drug substances: following the reference reading moved our
    agreement with FDA's published potency categories from 88% to 91%, which suggests
    phenolic beta-hydroxyls are what FDA actually counts in practice.
    """
    if ruleset.beta_hydroxyl == GUIDANCE:
        def matches(beta: Chem.Atom) -> bool:
            return (beta.GetHybridization() == Chem.HybridizationType.SP3
                    and _has_hydroxyl(beta))
    else:
        def matches(beta: Chem.Atom) -> bool:
            return (not _has_double_bond(beta)
                    and any(nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() > 0
                            for nb in beta.GetNeighbors()))

    return sum(1 for alpha in alphas
               if any(matches(b) for b in _beta_carbons(alpha, amine_n_idx)))


def _deactivating_factors(mol: Chem.Mol, site: SiteResult, amine_n_idx: int,
                          alphas: list[Chem.Atom], ruleset: Ruleset) -> list[Factor]:
    factors: list[Factor] = []

    if mol.HasSubstructMatch(rules.CARBOXYLIC_ACID):
        factors.append(Factor("Carboxylic acid group anywhere on molecule",
                              rules.CARBOXYLIC_ACID_SCORE, "Table B"))

    ring = _classify_ring(mol, amine_n_idx)
    if ring:
        score, label = ring
        factors.append(Factor(label, score, "Table B"))

    # Long chains on both sides, only for an acyclic N-nitroso group.
    if not mol.GetAtomWithIdx(amine_n_idx).IsInRing():
        if all(_chain_of_length_exists(mol, a, amine_n_idx,
                                       rules.CHAIN_MIN_LENGTH,
                                       rules.CHAIN_MAX_ATOMS_IN_ONE_RING)
               for a in alphas):
            factors.append(Factor(
                f"Chains of >={rules.CHAIN_MIN_LENGTH} consecutive non-hydrogen atoms "
                "on both sides of acyclic N-nitroso group",
                rules.LONG_CHAIN_BOTH_SIDES_SCORE, "Table B"))

    # Electron-withdrawing groups bonded to alpha-carbons. See cpca/ewg.py: the
    # guidance defers to a paper rather than listing the groups, so this reproduces
    # FDA's reference implementation, which enumerates ten patterns and counts
    # MATCHES rather than sides.
    matches = ewg.ewg_count(alphas, amine_n_idx)
    if matches:
        score = (rules.EWG_BOTH_SIDES_SCORE if matches > 1
                 else rules.EWG_ONE_SIDE_SCORE)
        side_word = "both sides" if matches > 1 else "one side"
        factors.append(Factor(
            f"Electron-withdrawing group bonded to alpha-carbon on {side_word}",
            score, "Table B",
            detail="matched: " + ", ".join(ewg.describe(alphas, amine_n_idx))))

    oh_sides = _beta_hydroxyl_sides(amine_n_idx, alphas, ruleset)
    if oh_sides:
        score = (rules.BETA_HYDROXYL_BOTH_SIDES_SCORE if oh_sides == 2
                 else rules.BETA_HYDROXYL_ONE_SIDE_SCORE)
        side_word = "both sides" if oh_sides == 2 else "one side"
        factors.append(Factor(
            f"Hydroxyl group bonded to beta-carbon on {side_word}", score, "Table B"))

    return factors


def _ewg_on_atom(mol: Chem.Mol, alpha: Chem.Atom) -> list[str]:
    """Electron-withdrawing groups bonded to this alpha-carbon. Table B footnote **
    excludes carboxylic acid, aryl and ketone, which are handled elsewhere."""
    neighbour_idxs = {nb.GetIdx() for nb in alpha.GetNeighbors()}
    acid_atoms = {i for match in mol.GetSubstructMatches(rules.CARBOXYLIC_ACID) for i in match}

    hits = []
    for name, patt in rules.EWG_PATTERNS:
        for match in mol.GetSubstructMatches(patt):
            if not (set(match) & neighbour_idxs):
                continue
            if set(match) & acid_atoms:
                continue  # carboxylic acid counted separately
            if _is_ketone(mol, match):
                continue  # ketone excluded: conflicting data
            hits.append(name)
            break
    return hits


def _is_ketone(mol: Chem.Mol, match: tuple[int, ...]) -> bool:
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        has_carbonyl_o = any(
            nb.GetAtomicNum() == 8
            and mol.GetBondBetweenAtoms(idx, nb.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
            for nb in atom.GetNeighbors())
        carbon_nbrs = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 6)
        if has_carbonyl_o and carbon_nbrs == 2:
            return True
    return False


def _activating_factors(mol: Chem.Mol, amine_n_idx: int,
                        alphas: list[Chem.Atom]) -> list[Factor]:
    factors: list[Factor] = []

    aryl = any(
        alpha.GetHybridization() == Chem.HybridizationType.SP3
        and any(nb.GetIsAromatic() for nb in alpha.GetNeighbors())
        for alpha in alphas)
    if aryl:
        factors.append(Factor("Aryl group bonded to alpha-carbon (benzylic substituent)",
                              rules.ARYL_ON_ALPHA_CARBON_SCORE, "Table C"))

    if _methyl_on_beta_carbon(amine_n_idx, alphas):
        factors.append(Factor("Methyl group bonded to beta-carbon",
                              rules.METHYL_ON_BETA_CARBON_SCORE, "Table C"))

    return factors


# Main entry point


def score_site(mol: Chem.Mol, site_id: int, amine_n_idx: int,
               nitroso_n_idx: int, ruleset: Ruleset | str | None = None) -> SiteResult:
    ruleset = policy.get(ruleset)
    site = SiteResult(site_id=site_id, amine_n_idx=amine_n_idx)
    alphas = _alpha_carbons(mol, amine_n_idx, nitroso_n_idx)

    counts = sorted(a.GetTotalNumHs() for a in alphas)
    site.alpha_h_counts = (counts[0], counts[1])

    # Figure 1 triage gates
    site.gate_any_alpha_h = sum(counts) > 0
    site.gate_multiple_alpha_h = max(counts) >= 2
    site.gate_tertiary_alpha_c = any(_is_tertiary_alpha_carbon(a, ruleset)
                                     for a in alphas)

    if not site.gate_any_alpha_h:
        site.triage_exit = "No hydrogens on either alpha-carbon (Figure 1, gate 1)"
        site.category = 5
    elif not site.gate_multiple_alpha_h:
        site.triage_exit = ("No alpha-carbon carries more than one hydrogen "
                            "(Figure 1, gate 2)")
        site.category = 5
    elif site.gate_tertiary_alpha_c:
        site.triage_exit = "Tertiary alpha-carbon present (Figure 1, gate 3)"
        site.category = 5

    if site.category == 5:
        site.ai_ng_day = rules.CATEGORY_AI_NG_DAY[5]
        return site

    # Table A
    # Divergence 1 of 4: which lookup table. See cpca/policy.py.
    table = (rules.ALPHA_HYDROGEN_SCORE if ruleset.alpha_h_table == GUIDANCE
             else rules.REFERENCE_ALPHA_HYDROGEN_SCORE)

    pair = (counts[0], counts[1])
    if pair not in table:
        site.notes.append(
            f"alpha-hydrogen pattern {pair[0]},{pair[1]} does not appear in Table A; "
            "this structure is outside the CPCA's defined scope and needs expert review")
        return site

    alpha_score = table[pair]
    detail = f"alpha-hydrogen counts {pair[0]},{pair[1]}"
    if pair not in rules.ALPHA_HYDROGEN_SCORE:
        detail += ("; this pattern is absent from the published Table A and is scored "
                   "here only because FDA's reference implementation scores it")

    # Divergence 2 of 4: the ethyl exception, Table A footnote $.
    if ruleset.ethyl_exception and pair == rules.ETHYL_EXCEPTION_PAIR:
        methylenes = [a for a in alphas if a.GetTotalNumHs() == 2]
        if any(_in_ethyl_group(a, amine_n_idx) for a in methylenes):
            alpha_score = rules.ETHYL_EXCEPTION_SCORE
            detail += "; methylene alpha-carbon is part of an ethyl group (footnote $)"

    site.alpha_h_score = alpha_score
    site.factors.append(Factor("alpha-Hydrogen Score", alpha_score, "Table A", detail))

    # Tables B and C
    site.factors.extend(
        _deactivating_factors(mol, site, amine_n_idx, alphas, ruleset))
    site.factors.extend(_activating_factors(mol, amine_n_idx, alphas))

    site.potency_score = sum(f.score for f in site.factors)
    site.category = rules.category_for_score(site.potency_score)
    site.ai_ng_day = rules.CATEGORY_AI_NG_DAY[site.category]
    return site


def score_smiles(smiles: str,
                 ruleset: Ruleset | str | None = None) -> CPCAResult:
    """Score one NDSRI. Returns a CPCAResult with a full factor breakdown.

    `ruleset` chooses what to follow where FDA's guidance and FDA's own code disagree:
    a Ruleset, a preset name ("balanced", "guidance", "reference"), or None for the
    default. See cpca/policy.py.
    """
    ruleset = policy.get(ruleset)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return CPCAResult(smiles=smiles, applicable=False,
                          ruleset=ruleset.name,
                          not_applicable_reason="SMILES could not be parsed")

    # Salts, hydrates and other counter-ions must not be scored. Table B's "carboxylic
    # acid group anywhere on molecule" is +3, so a maleate or citrate counter-ion would
    # otherwise deactivate a nitrosamine it is not even bonded to.
    #
    # Pick the largest fragment that actually CONTAINS an N-nitroso group, not simply
    # the largest fragment: a small nitrosamine paired with a bulky counter-ion (a
    # tosylate, a besylate, succinic acid) would otherwise have the counter-ion scored
    # in its place, which fails silently as "CPCA does not apply". Falling back to the
    # largest fragment keeps the old behaviour when no fragment carries the group.
    stripped: str | None = None
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(fragments) > 1:
        by_size = sorted(fragments, key=lambda m: m.GetNumHeavyAtoms(), reverse=True)
        with_core = [f for f in by_size if f.HasSubstructMatch(rules.NITROSAMINE_CORE)]
        parent = with_core[0] if with_core else by_size[0]

        others = [Chem.MolToSmiles(f) for f in fragments if f is not parent]
        note = ""
        if with_core and parent is not by_size[0]:
            note = (" The largest fragment was NOT the nitrosamine, so the "
                    "N-nitroso-bearing fragment was scored instead.")
        stripped = (f"Scored one fragment only; removed {len(others)} "
                    f"counter-ion/solvate fragment(s): "
                    f"{', '.join(sorted(set(others))[:4])}.{note}")
        mol = parent

    # Deuterium and other labels must not change the alpha-hydrogen count. See
    # _strip_isotopes: left in place they send labelled standards of potent
    # nitrosamines to Category 5.
    mol, isotope_note = _strip_isotopes(mol)

    matches = mol.GetSubstructMatches(rules.NITROSAMINE_CORE)
    if not matches:
        return CPCAResult(
            smiles=smiles, applicable=False, ruleset=ruleset.name,
            warnings=[n for n in (stripped, isotope_note) if n],
            not_applicable_reason=(
                "No N-nitrosamine group found (needs a carbon on both sides of an "
                "N-nitroso group). CPCA does not apply."))

    result = CPCAResult(smiles=smiles, applicable=True, ruleset=ruleset.name)
    if stripped:
        result.warnings.append(stripped)
    if isotope_note:
        result.warnings.append(isotope_note)

    # Deduplicate: the SMARTS matches once per ordering of the two alpha-carbons.
    seen: set[int] = set()
    sites: list[tuple[int, int]] = []
    for match in matches:
        amine_n, nitroso_n = match[0], match[3]
        if amine_n in seen:
            continue
        seen.add(amine_n)
        sites.append((amine_n, nitroso_n))

    for i, (amine_n, nitroso_n) in enumerate(sites, start=1):
        amine = mol.GetAtomWithIdx(amine_n)

        if amine.GetIsAromatic():
            result.warnings.append(
                f"Site {i}: N-nitroso group is within an aromatic ring; CPCA does not "
                "apply to this site (section IV.C)")
            continue

        alphas = _alpha_carbons(mol, amine_n, nitroso_n)
        excluded = [a for a in alphas if _double_bonded_to_heteroatom(mol, a)]
        if excluded:
            result.warnings.append(
                f"Site {i}: an alpha-carbon is double bonded to a heteroatom "
                "(N-nitrosamide / urea / guanidine); CPCA does not apply to this site "
                "(section IV.C)")
            continue

        result.sites.append(score_site(mol, i, amine_n, nitroso_n, ruleset))

    if len(sites) > 2:
        result.warnings.append(
            f"{len(sites)} N-nitroso groups found. For more than two, FDA directs the "
            "applicant to contact the Agency (section IV.C, footnote 26)")

    governing = result.governing_site
    if governing is None:
        result.applicable = False
        notes = [n for s in result.sites for n in s.notes]
        result.not_applicable_reason = (
            "; ".join(notes) if notes
            else "No N-nitroso site could be scored; see warnings")
        return result

    # Section IV.C: with more than one N-nitroso group, the one with the highest
    # predicted potency (lowest category number) governs the whole molecule.
    result.category = governing.category
    result.ai_ng_day = governing.ai_ng_day
    return result


def _strip_isotopes(mol: Chem.Mol) -> tuple[Chem.Mol, str | None]:
    """Score an isotope-labelled structure as its unlabelled form.

    RDKit represents a deuterium as an explicit hydrogen ATOM in the graph, not as an
    implicit hydrogen count, so a CD3 group reports GetTotalNumHs() == 0. Figure 1's
    first triage gate, "no hydrogens on either alpha-carbon", then fires and the
    compound is sent to Category 5, 1500 ng/day. For labelled standards of known potent
    nitrosamines that is wrong by up to 57-fold in the permissive direction:
    N-nitrosodipropylamine-D14 scored Category 5 where the unlabelled compound is
    Category 1.

    The CPCA has no provision for isotope labelling, and a CD3 group is chemically a
    methyl, so the labels are removed before scoring and the report says so. The
    kinetic isotope effect on alpha-hydroxylation is real but is not part of the CPCA
    and is not modelled here; a labelled internal standard is in any case not the thing
    a limit is being set for.
    """
    labelled = [a for a in mol.GetAtoms() if a.GetIsotope()]
    if not labelled:
        return mol, None

    kinds = sorted({f"{a.GetIsotope()}{a.GetSymbol()}" for a in labelled})
    plain = Chem.Mol(mol)
    for atom in plain.GetAtoms():
        atom.SetIsotope(0)
    try:
        plain = Chem.RemoveHs(plain)
        Chem.SanitizeMol(plain)
    except Exception:  # noqa: BLE001  - keep the labelled structure rather than fail
        return mol, ("Structure carries isotope labels that could not be removed; "
                     "alpha-hydrogen counts may be wrong.")

    note = (f"Isotope labels removed before scoring ({len(labelled)} atom(s): "
            f"{', '.join(kinds)}); the structure is scored as its unlabelled form. "
            "The CPCA has no provision for isotope labelling and does not model the "
            "kinetic isotope effect.")
    return plain, note


def _double_bonded_to_heteroatom(mol: Chem.Mol, atom: Chem.Atom) -> bool:
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() in (7, 8, 16):
            return True
    return False
