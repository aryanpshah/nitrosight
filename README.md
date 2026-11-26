# NitroSight

NitroSight works out how potent a nitrosamine impurity is likely to be, and
shows you how it reached that answer.

You give it a structure. It gives you back a potency category from 1 to 5, the
acceptable daily intake limit that goes with that category, and a list of every
factor that contributed to the score with the table of the FDA guidance each one
came from. You can give it a picture of a structure, a SMILES string, a compound
name, or a CAS number.

This is the software behind the paper "From chemical structure images to
explainable nitrosamine potency categorization".

## What problem this solves

Regulators need a limit for every nitrosamine impurity in a medicine. Most of
these impurities have never been tested in animals, so the limit is predicted
from the structure instead, using a method called the Carcinogenic Potency
Categorization Approach (CPCA).

Doing that by hand is slow. Existing software needs the structure already typed
out as a SMILES string, which is often not how you have it. And, as this project
found, the existing tools do not always agree with each other.

## How the scoring works

The method asks three questions first. Does either carbon next to the N-nitroso
group carry any hydrogens? Does either carry more than one? Is either a tertiary
carbon? If the answers rule out the metabolic step that makes nitrosamines
harmful, the compound goes straight to Category 5, the least restrictive.

Otherwise the compound gets a score. It starts from how many hydrogens sit on
each neighbouring carbon, then goes up for features that make activation harder
and down for features that make it easier. The total maps to a category, and the
category maps to a daily limit between 26.5 and 1500 nanograms.

NitroSight shows you every step of that sum rather than just the answer.

## Quick start

```
pip install -r requirements.txt
python score.py "CCN(N=O)CC"
```

That gives you:

```
RESULT: Potency Category 1    Recommended AI = 26.5 ng/day

N-nitroso site 1 (amine N at atom 2)
  alpha-hydrogen counts (lowest first): 2,2

  Factors:
    +1   alpha-Hydrogen Score  [Table A]
           alpha-hydrogen counts 2,2

  Potency Score = (+1) = 1
  -> Potency Category 1, AI = 26.5 ng/day
```

The `[Table A]` tag is the point. Every line tells you which part of the FDA
guidance it came from, so you can check the result against the document by hand.

If a molecule has more than one N-nitroso group, each one is scored separately
and the most potent one sets the limit. All of them appear in the report.

To see the window instead:

```
python run_gui.py
```

## What you need installed

Python 3.11 or newer, and the packages in `requirements.txt`.

| Package | What it is for |
|---|---|
| rdkit | Reading structures and scoring them. Nothing works without it. |
| PyQt6, Pillow | The desktop window |
| openpyxl | The Excel outputs in `tools/` |
| pymupdf, pypdf | Reading the guidance PDFs in `tools/` |

If you only want to score structures from the command line, `pip install rdkit`
is enough.

Reading structures from images needs DECIMER. It pulls in TensorFlow, so give it
its own environment rather than installing it alongside everything else:

```
py -3.13 -m venv .venv-osr
.venv-osr\Scripts\python -m pip install decimer
```

NitroSight finds that environment on its own. If it is missing, the image button
is greyed out and everything else still works.

## Other things you can run

```
python score.py "<SMILES>" --json              machine-readable output
python score.py "<SMILES>" --rules guidance    use a different rule reading
python score.py "<SMILES>" --compare-rules     all three readings side by side
```

## Which rules it follows

FDA published a guidance document and a reference program, and the two disagree
on four scoring rules. One of them is a rule the document requires that the
program computes and then never uses.

Because both are official, NitroSight treats the choice as a setting rather than
deciding for you. There are three options:

| Setting | What it does |
|---|---|
| `balanced` (default) | Follows the document where it is clear and the program contradicts it. Follows the program where the document is vague and the program settles it. |
| `guidance` | Follows the published document throughout. |
| `reference` | Reproduces FDA's program throughout, including where it departs from the document. |

Every report says which setting produced it.

The four rules they disagree on:

| Rule | The document says | The program does |
|---|---|---|
| Table A patterns | 6 rows | 10 rows, adding one that scores compounds the document puts out of scope |
| Ethyl exception | Apply it | Computes it, then never uses it |
| Tertiary alpha-carbon | Carbon with three carbon neighbours | Carbon with no hydrogens and four bonds |
| Beta-hydroxyl | Only on an sp3 carbon | Also counts phenols |

The default follows the document on the first two and the program on the last
two. The reasoning is that the document wins when it is clear and the program
contradicts it, and the program wins when the document is vague and the program
settles it.

This changes the answer for about 1 structure in 80. To see whether it changes
yours, run `--compare-rules`.

## How the code is laid out

| Path | What is in it |
|---|---|
| `cpca/rules.py` | Only the numbers and definitions copied from the guidance, kept apart from the algorithm so they can be checked against the PDF |
| `cpca/scorer.py` | The scoring itself |
| `cpca/ewg.py` | Electron-withdrawing group patterns, the least precisely specified part of the method |
| `cpca/policy.py` | The three rule settings described above |
| `cpca/report.py` | Turns a result into a readable report |
| `resolve/` | Looking up structures by name or CAS number, and building likely impurities from a parent drug |
| `osr.py` | Reading structures from images |
| `gui/` | The desktop window |
| `tools/` | Validation and analysis scripts, listed below |

Scoring is deterministic. The same structure always gives the same category, and
no machine learning is involved in it. Image recognition is the one learned
component, and its output is always shown to you before anything is scored.

## Reproducing the numbers in the paper

Each claim in the paper comes from one of these scripts.

| Claim in the paper | Command |
|---|---|
| FDA's four worked examples, 4 of 4 | `python test_fda_examples.py` |
| Agreement with FDA's reference tool | `python tools/benchmark_fda.py` |
| The same, under each rule setting | `python tools/benchmark_fda.py --all-rulesets` |
| Three tools compared | `python tools/benchmark_tools.py` |
| Agreement with FDA's published categories | `python tools/validate_vs_fda.py` |
| Agreement with manual scoring | `python tools/validate_vs_manual.py` |
| Checking for information leakage | `python tools/holdout.py` |
| Image recognition accuracy | `python tools/osr_accuracy.py` |
| Catalogue results and the priority list | `python tools/export_report.py` |
| Rescore without re-fetching from PubChem | `python tools/rescore.py` |

Each script writes its results into `data/`, so you can check any number in the
paper by rerunning the script and reading the file it produces.

These scripts overwrite the committed results in place. That is deliberate — the
file you get back is the one the paper cites — but it means a rerun replaces the
published artefact rather than writing a copy beside it. Work on a clean
checkout, or run `git diff data/` afterwards to see exactly what moved.
`rescore.py` takes `--dry-run` if you only want to know whether anything would
change.

One script is not reproducible by rerunning: `tools/fetch_structures.py` queries
PubChem live, so its output depends on what PubChem returns that day.
`tools/rescore.py` exists to re-score the committed structures without
refetching, and is the one to use when checking the scoring itself.

Tests: `test_cpca.py`, `test_fda_examples.py`, `test_policy.py`,
`test_pipeline.py`, `test_nitrosate.py`. All five run from a clean clone with no
setup beyond `pip install -r requirements.txt`, and all five should pass.

Some scripts in `tools/` need material that is not redistributed, for the
reasons below, and will not run from a clean clone:

| Script | Needs |
|---|---|
| `export_doubtful.py`, `count_doubtful.py` | the vendor workbook and the structure images |
| `extract_workbook_images.py` | the vendor workbook |
| `run_osr.py` | the structure images and their index, plus the DECIMER environment |
| `benchmark_fda.py` (full run) | FDA's Java tool built in `external/` |
| `benchmark_tools.py` | the Novartis tool cloned into `external/` |

They are included because they show how the data in `data/` was produced. What
they generated is committed, so `osr_accuracy.py`, `validate_vs_fda.py`,
`validate_vs_manual.py`, `holdout.py`, `export_report.py`, `rescore.py` and
`benchmark_fda.py --compare` all run from a clean clone against those cached
outputs.

The two comparison scripts need the other tools, which are not included here:

```
mkdir external && cd external
git clone https://github.com/FDA/featurize-nitrosamines
git clone https://github.com/Novartis/NDSRIs_in_silico_tool
```

FDA's tool is Java and needs Maven to build. Novartis's tool is Python. Note
that its `CPCA_rules.py` cannot be imported as published, because its licence
header has no comment markers on it; `tools/benchmark_tools.py` works around
that in memory and leaves their file alone.

## What is not in this repository, and why

**The source workbook and the structure images.** The compound list and the
drawings come from a commercial catalogue and belong to its publisher, so they
are not redistributed. Everything derived from them is here, including every
resolved structure, every score and every per-compound image recognition result,
so all the published numbers can still be checked.

**The guidance PDFs.** FDA documents are public domain, but the ICH and EMA ones
are not. Links are below.

## Limits on what this tool tells you

A potency category is a prediction from structure. It is not a measurement, and
it is not an estimate of risk to a patient.

A compound appearing in a supplier catalogue means a reference standard exists.
It does not mean the impurity has ever been found in a medicine.

Image recognition gets the category right for 38 of the 39 compounds whose
structure we could check against a verified source. That sounds better than it
is. Thirty-nine is a small sample, and across every compound we could check it
is right about 84 percent of the time, so roughly one in six would be wrong if
nobody checked. It does worst on large molecules and ones with several rings.
Always look at the structure it shows you.

This is a screening aid for a qualified assessor, not a regulatory
determination.

## Sources

- FDA, Recommended Acceptable Intake Limits for Nitrosamine Drug Substance-Related Impurities, August 2023: https://www.fda.gov/media/170794/download
- FDA, CDER nitrosamine impurity acceptable intake limits: https://www.fda.gov/regulatory-information/search-fda-guidance-documents/cder-nitrosamine-impurity-acceptable-intake-limits
- EMA, Appendix 2, Carcinogenic Potency Categorisation Approach: https://www.ema.europa.eu/en/documents/other/appendix-2-carcinogenic-potency-categorisation-approach-n-nitrosamines_en.pdf
- Kruhlak NL et al. Regul Toxicol Pharmacol 2024;150:105640. https://doi.org/10.1016/j.yrtph.2024.105640
- FDA Featurize-Nitrosamines: https://github.com/FDA/featurize-nitrosamines
- Novartis NDSRIs in silico tool: https://github.com/Novartis/NDSRIs_in_silico_tool
- DECIMER: https://decimer.ai

## Citing this

If you use NitroSight, please cite the paper. [CITATION.cff](CITATION.cff) has
the details in a form GitHub and most reference managers can read directly.

Each release is archived on Zenodo with its own DOI. Cite the version you
actually used; the concept DOI below always resolves to the latest one.

<!-- ZENODO_BADGE -->


## Licence

MIT. See [LICENSE](LICENSE).
