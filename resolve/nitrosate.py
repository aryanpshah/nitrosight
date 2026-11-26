"""
Build NDSRI structures from their parent drug.

73% of the catalogue has no CAS number and carries vendor names ("Abacavir Nitroso 1")
that no public database resolves. But an NDSRI is by definition the parent drug with an
N-nitroso group on one of its amines - so given the parent, the plausible impurities can
be enumerated directly.

Where a parent has exactly one nitrosatable site the product is unambiguous. Where it
has several, every candidate is generated and the caller is told the result is ambiguous;
for screening, the safest reading is the worst-case category across all candidates.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from itertools import combinations

from rdkit import Chem


# "Desmethyl" in a vendor name is not decoration: a tertiary amine cannot be
# nitrosated directly, so an impurity named desmethyl-X must have lost a methyl
# first. The absence of the word is the ordinary reading of an NDSRI name -
# direct nitrosation of the parent's own secondary amine.
_DESMETHYL_IN_NAME = re.compile(r"des-?methyl|de-?methylated", re.IGNORECASE)


def implied_route(impurity_name: str) -> str | None:
    """The formation route a vendor name implies, or None if there is no name.

    This is an inference from a vendor's wording, not a structural
    determination, so a caller that uses it must NOT raise its confidence on
    the strength of it. Its value is that it often picks the right structure
    out of several candidates, which is what makes the InChIKey and any
    comparison against a regulator's table refer to the right molecule.
    """
    if not impurity_name or not impurity_name.strip():
        return None
    if _DESMETHYL_IN_NAME.search(impurity_name):
        return "tertiary_demethylation"
    return "secondary_amine"


def name_implied_index(impurity_name: str, routes: list[str]) -> int | None:
    """Index of the ONE candidate whose route the vendor name implies.

    `routes` is the formation route of each candidate, in candidate order.
    Returns None when the name implies nothing, or when the implied route
    matches zero candidates or more than one: a name that does not single out
    one structure must not be treated as if it did.
    """
    route = implied_route(impurity_name)
    if route is None:
        return None
    hits = [i for i, r in enumerate(routes) if r == route]
    return hits[0] if len(hits) == 1 else None


@dataclass
class NitrosoCandidate:
    smiles: str
    site_atom_idxs: tuple[int, ...]
    site_description: str
    route: str = "secondary_amine"


def _is_nitrosatable(atom: Chem.Atom) -> bool:
    """A secondary amine that can form a stable N-nitrosamine.

    Excludes primary amines (give diazonium, not nitrosamines), tertiary amines (need
    nitrosative dealkylation - a different route), aromatic nitrogens, and amide /
    sulfonamide nitrogens (give N-nitrosamides, which CPCA explicitly excludes).
    """
    if atom.GetAtomicNum() != 7 or atom.GetIsAromatic():
        return False
    if atom.GetFormalCharge() != 0 or atom.GetTotalNumHs() != 1:
        return False

    carbons = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 6]
    if len(carbons) != 2 or len(atom.GetNeighbors()) != 2:
        return False

    # Reject amide / sulfonamide / related: neighbour double-bonded to a heteroatom.
    for nb in atom.GetNeighbors():
        if nb.GetAtomicNum() not in (6, 16):
            continue
        for bond in nb.GetBonds():
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
            if bond.GetOtherAtom(nb).GetAtomicNum() in (7, 8, 16):
                return False
    return True


def find_nitrosatable_sites(mol: Chem.Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if _is_nitrosatable(a)]


def _is_methyl(atom: Chem.Atom) -> bool:
    return (atom.GetAtomicNum() == 6
            and atom.GetDegree() == 1
            and atom.GetTotalNumHs() == 3)


def find_demethylation_sites(mol: Chem.Mol) -> list[tuple[int, int]]:
    """Tertiary amines bearing an N-methyl, as (amine N index, methyl C index).

    Tertiary amines do not nitrosate directly - they undergo nitrosative dealkylation,
    losing an alkyl group first. Loss of a methyl is the common case and gives the
    "N-nitroso-desmethyl-X" impurities that appear throughout the catalogue.
    """
    sites: list[tuple[int, int]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetIsAromatic():
            continue
        if atom.GetFormalCharge() != 0 or atom.GetTotalNumHs() != 0:
            continue
        neighbours = atom.GetNeighbors()
        carbons = [nb for nb in neighbours if nb.GetAtomicNum() == 6]
        if len(carbons) != 3 or len(neighbours) != 3:
            continue

        # Same amide / sulfonamide exclusion as for secondary amines.
        conjugated = False
        for nb in neighbours:
            if nb.GetAtomicNum() not in (6, 16):
                continue
            for bond in nb.GetBonds():
                if (bond.GetBondType() == Chem.BondType.DOUBLE
                        and bond.GetOtherAtom(nb).GetAtomicNum() in (7, 8, 16)):
                    conjugated = True
        if conjugated:
            continue

        for methyl in (c for c in carbons if _is_methyl(c)):
            sites.append((atom.GetIdx(), methyl.GetIdx()))
    return sites


def _describe_site(mol: Chem.Mol, idx: int) -> str:
    atom = mol.GetAtomWithIdx(idx)
    env = []
    for nb in atom.GetNeighbors():
        if nb.GetIsAromatic():
            env.append("aromatic C")
        elif nb.IsInRing():
            env.append("ring C")
        else:
            env.append("chain C")
    ring = " (in ring)" if atom.IsInRing() else ""
    return f"N{idx} between {' and '.join(env)}{ring}"


def _nitrosate(mol: Chem.Mol, ops: tuple[tuple[str, int, int | None], ...]) -> str | None:
    """Apply nitrosation operations. Atoms are added first (which leaves existing
    indices untouched), then any displaced methyls are removed highest-index first."""
    rw = Chem.RWMol(mol)
    to_remove: list[int] = []

    for _kind, n_idx, methyl_idx in ops:
        n = rw.GetAtomWithIdx(n_idx)
        nitroso_n = rw.AddAtom(Chem.Atom(7))
        oxygen = rw.AddAtom(Chem.Atom(8))
        rw.AddBond(n_idx, nitroso_n, Chem.BondType.SINGLE)
        rw.AddBond(nitroso_n, oxygen, Chem.BondType.DOUBLE)
        n.SetNumExplicitHs(0)
        n.SetNoImplicit(True)
        if methyl_idx is not None:
            to_remove.append(methyl_idx)

    for idx in sorted(to_remove, reverse=True):
        rw.RemoveAtom(idx)

    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:  # noqa: BLE001
        return None
    return Chem.MolToSmiles(out)


MAX_CANDIDATES = 12


def enumerate_ndsris(parent_smiles: str, n_groups: int = 1) -> list[NitrosoCandidate]:
    """All plausible ways of putting `n_groups` N-nitroso groups on the parent drug.

    Two formation routes are considered: direct nitrosation of a secondary amine, and
    nitrosative demethylation of a tertiary N-methyl amine. Candidates using only the
    direct route are returned first, being the more straightforward chemistry.
    """
    mol = Chem.MolFromSmiles(parent_smiles)
    if mol is None:
        return []

    ops: list[tuple[str, int, int | None]] = [
        ("secondary_amine", idx, None) for idx in find_nitrosatable_sites(mol)
    ] + [
        ("tertiary_demethylation", n_idx, me_idx)
        for n_idx, me_idx in find_demethylation_sites(mol)
    ]
    if len(ops) < n_groups:
        return []

    out: list[NitrosoCandidate] = []
    seen: set[str] = set()

    for combo in combinations(ops, n_groups):
        if len({op[1] for op in combo}) != n_groups:
            continue  # two operations on the same nitrogen
        smiles = _nitrosate(mol, combo)
        if not smiles or smiles in seen:
            continue
        seen.add(smiles)

        kinds = {op[0] for op in combo}
        out.append(NitrosoCandidate(
            smiles=smiles,
            site_atom_idxs=tuple(op[1] for op in combo),
            site_description="; ".join(_describe_site(mol, op[1]) for op in combo),
            route=kinds.pop() if len(kinds) == 1 else "mixed",
        ))

    # Direct nitrosation first, then demethylation, then mixed.
    order = {"secondary_amine": 0, "tertiary_demethylation": 1, "mixed": 2}
    out.sort(key=lambda c: order.get(c.route, 3))
    return out[:MAX_CANDIDATES]


def groups_from_name(impurity_name: str) -> int:
    """How many N-nitroso groups does the vendor's name imply?"""
    low = impurity_name.lower()
    if "trinitroso" in low:
        return 3
    if "dinitroso" in low or "di-nitroso" in low or "bis-nitroso" in low:
        return 2
    return 1
