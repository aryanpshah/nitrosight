"""
Draw the impurity, highlighting the atoms the score actually depends on.

The N-nitroso group and its alpha- and beta-carbons are what CPCA counts, so colouring
them turns the factor list into something a chemist can check at a glance.
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

NITROSO = Chem.MolFromSmarts("[NX3](-[#6])(-[#6])-[NX2]=[OX1]")

COLOUR_NITROSO = (0.93, 0.33, 0.31)   # red    - the N-nitroso group itself
COLOUR_ALPHA = (0.96, 0.66, 0.26)     # amber  - alpha-carbons (hydrogen count)
COLOUR_BETA = (0.40, 0.68, 0.90)      # blue   - beta-carbons (hydroxyl, methyl)


def _highlight_map(mol: Chem.Mol) -> tuple[list[int], dict[int, tuple]]:
    atoms: dict[int, tuple] = {}

    for match in mol.GetSubstructMatches(NITROSO):
        amine_n, a1, a2, nitroso_n = match[0], match[1], match[2], match[3]
        atoms[amine_n] = COLOUR_NITROSO
        atoms[nitroso_n] = COLOUR_NITROSO
        for nb in mol.GetAtomWithIdx(nitroso_n).GetNeighbors():
            if nb.GetAtomicNum() == 8:
                atoms[nb.GetIdx()] = COLOUR_NITROSO

        for alpha in (a1, a2):
            atoms.setdefault(alpha, COLOUR_ALPHA)
            for nb in mol.GetAtomWithIdx(alpha).GetNeighbors():
                if nb.GetAtomicNum() == 6 and nb.GetIdx() != amine_n:
                    atoms.setdefault(nb.GetIdx(), COLOUR_BETA)

    return list(atoms), atoms


def render_png(smiles: str, width: int = 460, height: int = 330,
               highlight: bool = True) -> bytes | None:
    """Return PNG bytes for this structure, or None if it cannot be drawn."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Score and draw the parent fragment, consistent with the scorer.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

    Chem.rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    opts.addStereoAnnotation = False
    opts.bondLineWidth = 2

    if highlight:
        ids, colours = _highlight_map(mol)
        drawer.DrawMolecule(mol, highlightAtoms=ids, highlightAtomColors=colours)
    else:
        drawer.DrawMolecule(mol)

    drawer.FinishDrawing()
    return drawer.GetDrawingText()


LEGEND = [
    ("N-nitroso group", COLOUR_NITROSO),
    ("α-carbons (hydrogen count)", COLOUR_ALPHA),
    ("β-carbons (hydroxyl / methyl)", COLOUR_BETA),
]


def rgb_css(colour: tuple[float, float, float]) -> str:
    r, g, b = (int(round(c * 255)) for c in colour)
    return f"rgb({r},{g},{b})"
