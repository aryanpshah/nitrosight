"""
Electron-withdrawing groups on the alpha-carbon (Table B).

This is the least precisely specified rule in the CPCA guidance. Table B's footnote **
says only that EWGs exclude carboxylic acid, aryl and ketone, and defers to
Cross KP & Ponting DJ (2021), Comput Toxicol 20:100186, rather than listing them.

FDA's reference implementation settles it. featurize-nitrosamines
(FeaturizeNitrosamine.calculateFeatures) sums ten separate, individually enumerated
patterns into one count, and then:

    exactly one match  -> "EWG on one side",   +1
    more than one      -> "EWG on both sides", +2

Note that FDA counts MATCHES, not sides: two EWGs on a single alpha-carbon score as
"both sides". That is reproduced here deliberately, so the two implementations agree.

Also note FDA computes an eleventh pattern, ewg_aromaticBeta (N-CH2-[O,N,S]-aromatic),
but does NOT add it to the total. That omission is reproduced too.

Every predicate below mirrors the corresponding method in FDA's AtomTest class; the
Java name is given in each docstring so the two can be read side by side.
"""

from __future__ import annotations

from rdkit import Chem

HALOGENS = {9, 17, 35, 53}       # F, Cl, Br, I


# Atom predicates, mirroring gov.fda.gsrs.chem.util.AtomTest


def _h(atom: Chem.Atom) -> int:
    """AtomTest.getHCount"""
    return atom.GetTotalNumHs()


def _carbons(atom: Chem.Atom) -> int:
    """AtomTest.countCarbon"""
    return sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 6)


def _double_bonded(atom: Chem.Atom) -> list[Chem.Atom]:
    return [b.GetOtherAtom(atom) for b in atom.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE]


def _has_double_bond(atom: Chem.Atom) -> bool:
    return bool(_double_bonded(atom))


def is_double_bonded_oxygen(atom: Chem.Atom) -> bool:
    """AtomTest.isDoubleBondedOxygen"""
    return atom.GetAtomicNum() == 8 and _has_double_bond(atom)


def is_carbonyl(atom: Chem.Atom) -> bool:
    """AtomTest.isCarbonyl"""
    return atom.GetAtomicNum() == 6 and any(
        is_double_bonded_oxygen(nb) for nb in atom.GetNeighbors())


def is_sulfone(atom: Chem.Atom) -> bool:
    """AtomTest.isSulfone: sulfur with exactly two doubly bonded oxygens."""
    return atom.GetAtomicNum() == 16 and sum(
        1 for nb in atom.GetNeighbors() if is_double_bonded_oxygen(nb)) == 2


def _has_oh(atom: Chem.Atom) -> bool:
    return any(nb.GetAtomicNum() == 8 and _h(nb) >= 1 and nb.GetDegree() == 1
               for nb in atom.GetNeighbors())


def is_nitrate(atom: Chem.Atom) -> bool:
    """AtomTest.isNitrate: nitrogen, no hydrogens, one doubly bonded oxygen, has an OH."""
    return (atom.GetAtomicNum() == 7 and _h(atom) == 0
            and sum(1 for nb in atom.GetNeighbors()
                    if is_double_bonded_oxygen(nb)) == 1
            and _has_oh(atom))


def is_quat_amine(atom: Chem.Atom) -> bool:
    """AtomTest.isQuatAmine: nitrogen, no hydrogens, four bonds."""
    return (atom.GetAtomicNum() == 7 and _h(atom) == 0
            and atom.GetDegree() == 4)


def is_cyanide(atom: Chem.Atom) -> bool:
    """AtomTest.isCyanide: the carbon of a nitrile."""
    if atom.GetAtomicNum() != 6:
        return False
    return any(nb.GetAtomicNum() == 7
               and atom.GetOwningMol().GetBondBetweenAtoms(
                   atom.GetIdx(), nb.GetIdx()).GetBondType() == Chem.BondType.TRIPLE
               for nb in atom.GetNeighbors())


def is_double_bond_to_s_or_n(atom: Chem.Atom) -> bool:
    """AtomTest.isDoubleBondToSorN: carbon doubly bonded to nitrogen or sulfur."""
    if atom.GetAtomicNum() != 6:
        return False
    return any(nb.GetAtomicNum() in (7, 16) for nb in _double_bonded(atom))


def is_trifluoromethyl(atom: Chem.Atom) -> bool:
    return atom.GetAtomicNum() == 6 and sum(
        1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 9) == 3


# The ten patterns


def _ch1_or_ch2(atom: Chem.Atom) -> bool:
    """FDA gates several patterns on the alpha-carbon carrying one or two hydrogens."""
    return _h(atom) in (1, 2)


def _beta_carbonyl_derivative(alpha: Chem.Atom, amine_idx: int) -> int:
    """N-[CH2]-C(=[O,N,S])-([H],[C],[F,Cl,Br,I],[OC],[NH2],[NHC],[NC2]), ketone excluded."""
    if not _ch1_or_ch2(alpha):
        return 0

    count = 0
    for beta in alpha.GetNeighbors():
        if beta.GetAtomicNum() != 6 or beta.GetIdx() == amine_idx:
            continue
        if not any(nb.GetAtomicNum() in (7, 8, 16) for nb in _double_bonded(beta)):
            continue

        if _h(beta) == 1:                       # aldehyde-like
            count += 1
            continue

        mol = beta.GetOwningMol()
        for bond in beta.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            other = bond.GetOtherAtom(beta)
            sym = other.GetAtomicNum()
            ok = (
                sym == 1
                # a plain carbon substituent counts only when beta is NOT a carbonyl,
                # which is how FDA excludes ketones
                or (sym == 6 and _carbons(beta) == 2 and not is_carbonyl(beta))
                or sym in HALOGENS
                or (sym == 8 and _carbons(other) == 2)                  # ester
                or (sym == 7 and _h(other) == 2)                        # primary amide
                or (sym == 7 and _h(other) == 1 and _carbons(other) == 2)
                or (sym == 7 and _h(other) == 0 and _carbons(other) == 3)
            )
            if ok:
                count += 1
                break
        else:
            # beta carries implicit hydrogens only; FDA's bond scan sees those as
            # single bonds to hydrogen
            if _h(beta) > 0:
                count += 1
    return count


def _alpha_halogen(alpha: Chem.Atom) -> int:
    """N-[CH2]-[X]"""
    if not _ch1_or_ch2(alpha):
        return 0
    return sum(1 for nb in alpha.GetNeighbors() if nb.GetAtomicNum() in HALOGENS)


def _conjugated_gammas(alpha: Chem.Atom, amine_idx: int) -> list[Chem.Atom]:
    """Shared opening of the three conjugation patterns:
    N-[CH2]-[CH|C][=|#][CH|C]-  , returning the far end of the multiple bond."""
    if not _ch1_or_ch2(alpha):
        return []

    out = []
    for beta in alpha.GetNeighbors():
        if beta.GetAtomicNum() != 6 or beta.GetIdx() == amine_idx:
            continue
        if _h(beta) >= 2:
            continue
        if sum(1 for nb in beta.GetNeighbors() if nb.GetAtomicNum() != 1) != 2:
            continue
        for bond in beta.GetBonds():
            if bond.GetBondType() not in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                continue
            gamma = bond.GetOtherAtom(beta)
            if sum(1 for nb in gamma.GetNeighbors() if nb.GetAtomicNum() != 1) == 2:
                out.append(gamma)
    return out


def _aromatic_conjugated(alpha: Chem.Atom, amine_idx: int) -> int:
    """...-c : conjugated through to an aromatic carbon."""
    return sum(1 for gamma in _conjugated_gammas(alpha, amine_idx)
               if any(nb.GetIsAromatic() and nb.GetAtomicNum() == 6
                      for nb in gamma.GetNeighbors()))


def _michael_acceptor_imine(alpha: Chem.Atom, amine_idx: int) -> int:
    return sum(
        1 for gamma in _conjugated_gammas(alpha, amine_idx)
        if any(any(is_double_bond_to_s_or_n(far) or is_double_bonded_oxygen(far)
                   for far in nb.GetNeighbors())
               for nb in gamma.GetNeighbors()))


def _michael_acceptor_cyanide(alpha: Chem.Atom, amine_idx: int) -> int:
    return sum(1 for gamma in _conjugated_gammas(alpha, amine_idx)
               if any(is_cyanide(nb) for nb in gamma.GetNeighbors()))


def _neighbour_matches(alpha: Chem.Atom, predicate) -> int:
    return sum(1 for nb in alpha.GetNeighbors() if predicate(nb))


def ewg_count(alphas: list[Chem.Atom], amine_idx: int) -> int:
    """Total EWG matches around this N-nitroso group, as FDA sums them."""
    total = 0
    for alpha in alphas:
        total += _beta_carbonyl_derivative(alpha, amine_idx)
        total += _alpha_halogen(alpha)
        total += _aromatic_conjugated(alpha, amine_idx)
        total += _michael_acceptor_imine(alpha, amine_idx)
        total += _michael_acceptor_cyanide(alpha, amine_idx)
        total += _neighbour_matches(alpha, is_sulfone)
        total += _neighbour_matches(alpha, is_nitrate)
        total += _neighbour_matches(alpha, is_quat_amine)
        total += _neighbour_matches(alpha, is_cyanide)
        total += _neighbour_matches(
            alpha, lambda nb: nb.GetAtomicNum() == 6 and any(
                is_trifluoromethyl(n2) for n2 in nb.GetNeighbors()))
    return total


def describe(alphas: list[Chem.Atom], amine_idx: int) -> list[str]:
    """Which patterns fired, for the report."""
    found: list[str] = []
    for alpha in alphas:
        checks = [
            ("carbonyl derivative on beta-carbon",
             _beta_carbonyl_derivative(alpha, amine_idx)),
            ("halogen on alpha-carbon", _alpha_halogen(alpha)),
            ("conjugated to an aromatic ring",
             _aromatic_conjugated(alpha, amine_idx)),
            ("Michael acceptor (imine/carbonyl)",
             _michael_acceptor_imine(alpha, amine_idx)),
            ("Michael acceptor (nitrile)",
             _michael_acceptor_cyanide(alpha, amine_idx)),
            ("sulfone", _neighbour_matches(alpha, is_sulfone)),
            ("nitrate", _neighbour_matches(alpha, is_nitrate)),
            ("quaternary nitrogen", _neighbour_matches(alpha, is_quat_amine)),
            ("nitrile", _neighbour_matches(alpha, is_cyanide)),
            ("trifluoromethyl", _neighbour_matches(
                alpha, lambda nb: nb.GetAtomicNum() == 6 and any(
                    is_trifluoromethyl(n2) for n2 in nb.GetNeighbors()))),
        ]
        found += [name for name, n in checks if n]
    return sorted(set(found))
