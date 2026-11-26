# NitroSight

NitroSight is an open-source tool that predicts how potent a nitrosamine impurity is likely to be. It uses the FDA's Carcinogenic Potency Categorization Approach (CPCA).

You can give it a picture of a structure, a SMILES string, a compound name, or a CAS number. It gives back a potency category, the acceptable intake (AI) limit for that category, and every step of the score so you can check it.

This repository goes with the paper *NitroSight: from chemical structure images to explainable nitrosamine potency categorization*.

## Features

- CPCA scoring and potency categories
- Reading structures from pictures with DECIMER
- Looking up compounds by name or CAS number
- Building the likely nitrosamines of a parent drug
- Every score points to the FDA guidance table it came from
- Molecules with more than one N-nitroso group
- A desktop window and a command line
- Three settings for the rules where FDA's two sources disagree

## Installation

You need Python 3.11 or newer.

```bash
pip install -r requirements.txt
```

| Package | What it is for |
|---|---|
| RDKit | Reading structures and scoring them |
| PyQt6 | The desktop window |
| openpyxl | The Excel files written by the scripts in `tools/` |

If you only want to score from the command line, this is enough:

```bash
pip install rdkit
```

### Reading structures from pictures (optional)

This needs DECIMER. DECIMER brings in TensorFlow, so give it its own Python environment:

```bash
py -3.13 -m venv .venv-osr
.venv-osr\Scripts\python -m pip install decimer
```

NitroSight finds this environment on its own. Without it, everything else still works.

## Usage

### Desktop window

```bash
python run_gui.py
```

### Command line

Score a compound from its SMILES string:

```bash
python score.py "CCN(N=O)CC"
```

Part of the output:

```text
RESULT: Potency Category 1    Recommended AI = 26.5 ng/day

N-nitroso site 1 (amine N at atom 2)
  alpha-hydrogen counts (lowest first): 2,2

  Factors:
    +1   alpha-Hydrogen Score  [Table A]
           alpha-hydrogen counts 2,2

  Potency Score = (+1) = 1
  -> Potency Category 1, AI = 26.5 ng/day
```

The `[Table A]` tag tells you which part of the FDA guidance each line came from.

Other options:

```bash
python score.py "<SMILES>" --json             # machine-readable output
python score.py "<SMILES>" --rules guidance   # use a different rule setting
python score.py "<SMILES>" --compare-rules    # all three settings side by side
```

## How the scoring works

The CPCA puts a nitrosamine in one of five categories. The limits go from 26.5 ng/day (Category 1, most potent) to 1500 ng/day (Categories 4 and 5).

The score is built from:

- how many hydrogens sit on the carbons next to the N-nitroso group
- features that make the compound more potent
- features that make it less potent
- rules that say when the method does not apply

If a molecule has more than one N-nitroso group, each one is scored and the most potent one sets the limit.

The scoring has no machine learning in it. The same structure always gives the same answer. Reading a picture is the only learned step, and NitroSight always shows you the structure it read before scoring it.

## Rule settings

FDA published two things that define the CPCA: a guidance document and a reference program. They disagree on four rules. NitroSight lets you choose which one to follow.

| Setting | What it does |
|---|---|
| `balanced` (default) | Follows the document where it is clear, and the program where the document is vague |
| `guidance` | Follows the FDA guidance document |
| `reference` | Follows the FDA reference program |

The four rules they disagree on:

| Rule | FDA guidance | FDA reference program |
|---|---|---|
| Table A patterns | 6 patterns | 10 patterns |
| Ethyl exception | Applied | Worked out but never used |
| Tertiary alpha-carbon | Three carbon neighbours | No hydrogens and four bonds |
| Beta-hydroxyl | sp3 carbons only | Also counts phenols |

`balanced` follows the guidance on the first two rules and the program on the last two. Every report says which setting was used. To see if the setting changes your answer, use `--compare-rules`.

[docs/reference-implementation-divergences.md](docs/reference-implementation-divergences.md) explains each rule in more detail.

## Repository layout

```text
nitrosight/
├── cpca/
│   ├── rules.py       # numbers and definitions copied from the guidance
│   ├── scorer.py      # the scoring
│   ├── ewg.py         # electron-withdrawing groups
│   ├── policy.py      # the three rule settings
│   └── report.py      # the written report
├── resolve/           # name and CAS lookup, building nitrosamines from a parent drug
├── gui/               # the desktop window
├── tools/             # scripts that produce the results in the paper
├── data/              # the data and results
├── docs/              # where the guidance and the program disagree
├── osr.py             # reading structures from pictures
├── pipeline.py        # from any input to a final answer
├── score.py           # command line
├── run_gui.py         # starts the desktop window
└── test_*.py          # tests
```

## Reproducing the results

These scripts produce the results in the paper.

| Result | Command |
|---|---|
| FDA's four worked examples | `python test_fda_examples.py` |
| Agreement with FDA's reference program | `python tools/benchmark_fda.py --compare` |
| The same, under each rule setting | `python tools/benchmark_fda.py --all-rulesets` |
| Comparison of three tools | `python tools/benchmark_tools.py` |
| Agreement with FDA's published categories | `python tools/validate_vs_fda.py` |
| Agreement with manual scoring | `python tools/validate_vs_manual.py` |
| Hold-out check | `python tools/holdout.py` |
| Picture reading accuracy | `python tools/osr_accuracy.py` |
| Catalogue results | `python tools/export_report.py` |
| Rescoring the saved structures | `python tools/rescore.py` |

The results are saved in `data/`. Running a script overwrites its file there. Use a clean copy of the repository, or check what changed with:

```bash
git diff data/
```

To see what rescoring would change without saving anything:

```bash
python tools/rescore.py --dry-run
```

`tools/fetch_structures.py` downloads structures from PubChem, so its results can change over time. Use `tools/rescore.py` to rescore the saved structures instead.

### Tests

```bash
python test_cpca.py
python test_fda_examples.py
python test_policy.py
python test_pipeline.py
python test_nitrosate.py
```

All five run from a fresh copy after `pip install -r requirements.txt`.

### Other tools used in the comparison

The full comparison needs FDA's and Novartis's tools:

```bash
mkdir external
cd external
git clone https://github.com/FDA/featurize-nitrosamines
git clone https://github.com/Novartis/NDSRIs_in_silico_tool
```

FDA's tool needs Java and Maven. It is used by `python tools/benchmark_fda.py` with no options, which runs it again from scratch. Novartis's tool is in Python. Its `CPCA_rules.py` cannot be imported as published, so `tools/benchmark_tools.py` works around this without changing their file.

## Data availability

The compound list and structure pictures come from a commercial catalogue, so they are not included. Everything made from them is in `data/`: the structures, the scores, and the picture reading results. That means most scripts run without the catalogue.

These scripts need the original files:

| Script | Needs |
|---|---|
| `tools/extract_workbook_images.py` | the catalogue workbook |
| `tools/run_osr.py` | the structure pictures, their index, and DECIMER |

The guidance PDFs are not included. The FDA guidance is at https://www.fda.gov/media/170794/download.

## Limitations

- A CPCA category is a prediction from structure. It is not a measured potency.
- Picture reading can get structures wrong, more often for large molecules or ones with many rings. Always check the structure it shows you.
- A compound being in a supplier catalogue does not mean it has been found in a medicine.
- NitroSight is a screening tool for research. It does not replace an expert or a regulatory decision.

## Citation

If you use NitroSight, please cite the paper. [CITATION.cff](CITATION.cff) has the details.

The code is archived on Zenodo at https://doi.org/10.5281/zenodo.24572986. This link always points to the newest version. Version 1.0.0, the one described in the paper, is at https://doi.org/10.5281/zenodo.24572987.

## License

MIT. See [LICENSE](LICENSE).
