# Where the guidance and the reference implementation disagree

FDA published two things that define the CPCA: a guidance document, and a reference
program. On four scoring rules they do not say the same thing. Both are official, so
NitroSight does not pick a winner and hide the choice. The reading is a setting, every
report says which setting produced it, and this note records what each setting does and
why.

Sources:

- **Guidance** — FDA, *Recommended Acceptable Intake Limits for Nitrosamine Drug
  Substance-Related Impurities (NDSRIs)*, August 2023.
  <https://www.fda.gov/media/170794/download>. Figure 1 is on p.10, Figure 2 on p.11,
  and Appendix A (Tables A, B and C) on pp.18-20.
- **Reference code** — FDA, Featurize-Nitrosamines.
  <https://github.com/FDA/featurize-nitrosamines>. The relevant logic is in
  `FeaturizeNitrosamine.java` and `AtomTest.java`.
- EMA's Appendix 2 (EMA/451665/2023) carries the same tables as the FDA guidance and is
  used below where it corroborates a reading.

The three settings are defined in `cpca/policy.py`. `test_policy.py` is the executable
form of this document: every compound named here is a test case there.

| Rule | Guidance | Reference code | `guidance` | `reference` | `balanced` (default) |
|---|---|---|---|---|---|
| 1. Table A entries | 6 rows | 10 rows | guidance | code | **guidance** |
| 2. Ethyl exception | required | computed, never used | guidance | code | **guidance** |
| 3. Tertiary alpha-carbon | sp3, 3 carbon neighbours | 0 H and 4 bonds | guidance | code | **code** |
| 4. Beta-hydroxyl | sp3 beta-carbon only | phenols also count | guidance | code | **code** |

---

## 1. Table A has six rows, and the code uses ten

**Guidance.** Table A lists six alpha-hydrogen patterns, written as exact counts with
the lower count first: `0,2 -> 3`, `0,3 -> 2`, `1,2 -> 3`, `1,3 -> 3`, `2,2 -> 1`,
`2,3 -> 1`. Appendix A says an NDSRI "is expected to match exactly one of the
alpha-hydrogen definitions in Table A".

**Reference code.** `FeaturizeNitrosamine.java` builds a ten-entry lookup, adding
`0,0 -> 5`, `0,1 -> 4`, `1,1 -> 4` and `3,3 -> 1`.

**What follows.** Three of the four extra rows are unreachable: `0,0` fails Figure 1's
first triage gate, and `0,1` and `1,1` fail the second, so all three exit to Category 5
before Table A is consulted. Only `3,3` is reachable, and it is the one that matters: it
is the pattern of two methyl groups, so it is **N-nitrosodimethylamine (NDMA)** and its
relatives.

Read literally, the guidance therefore has no row for NDMA, and the `guidance` and
`balanced` settings decline to score it rather than guess:

```
$ python score.py "CN(C)N=O"
CPCA does not apply ...
$ python score.py "CN(C)N=O" --rules reference
RESULT: Potency Category 1    Recommended AI = 26.5 ng/day
```

This is deliberate. Declining is the honest answer when a structure falls outside the
published table, and silently inventing a row is how two implementations drift apart.
It is also worth stating plainly that this is a gap in Table A as published: the
archetypal nitrosamine, the one found in valsartan in 2018, does not match any row of
the table that is supposed to cover every NDSRI. NDMA has a compound-specific FDA limit
in practice, so nothing turns on it operationally, but any implementation has to decide
what to do here and the guidance does not say.

`balanced` follows the guidance because the document is explicit about matching exactly
one listed definition.

## 2. The ethyl exception is required, computed, and then discarded

**Guidance.** Table A's footnote `$`, attached to the `0,2` row, reads in full:

> A score of 3 applies when the methylene alpha-carbon is not part of an ethyl group.
> If the methylene alpha-carbon is part of an ethyl group, a score of 2 should be
> applied.

**Reference code.** The lookup maps `0,2` to 3 unconditionally. The program does compute
the ethyl condition, into a variable named `hasEthane`, and never reads the result; the
alpha-hydrogen score has already been fixed earlier in the same function.

**What follows.** For N-nitroso-N-ethylaniline the guidance gives Category 2 and a limit
of 100 ng/day; the reference code gives Category 3 and 400 ng/day. A four-fold
difference in an enforceable limit.

Two things indicate the footnote is intended rather than an artefact of the text. The
structure drawn in the `0,2` row of Table A is the propyl compound, not the ethyl one,
which is the natural way to illustrate a score of 3 without invoking the exception. And
EMA's Appendix 2 carries the identical footnote, so it is not a transcription slip
confined to one document.

`balanced` follows the guidance: the document is unambiguous and the code contradicts
it, and the dead variable in the code reads as an incomplete implementation rather than
a deliberate decision.

This finding has been reported to FDA as a public issue on the reference
implementation's repository:
<https://github.com/FDA/featurize-nitrosamines/issues/3>.

## 3. "Tertiary alpha-carbon" is defined twice

**Guidance.** Figure 1's footnote `*` defines the term:

> A tertiary alpha-carbon is defined as an alpha-carbon atom in an sp3 hybridization
> state, bonded to three other carbon atoms.

**Reference code.** `AtomTest` asks whether the atom has no hydrogens and four bonds.

**What follows.** The code's test is strictly broader. An alpha-carbon bonded to the
N-nitroso nitrogen, two carbons and one oxygen has no hydrogens and four bonds, so the
code calls it tertiary, while the guidance does not, because only two of its neighbours
are carbons. Triage gate 3 then sends the compound to Category 5 under the code's
reading and on to scoring under the guidance's. For `CC1CN(C(O1)(C)C)N=O` the guidance
gives Category 4 and the code gives Category 5.

**`balanced` follows the code here, and the reason is not that the guidance is unclear.**
It is not: the footnote is an explicit definition. The default follows the code because
FDA's own (Q)SAR Model Reporting Format names Featurize-Nitrosamines as the CPCA's
automated implementation, which makes it the de facto standard a submission will be
checked against. Users who want conformance to the document as written should select
`--rules guidance`, which is provided for exactly that purpose.

Note that the code's reading is the more permissive of the two: it sends more compounds
to Category 5 and its 1500 ng/day limit.

## 4. Beta-hydroxyl: whether a phenol counts

**Guidance.** Table B's footnote `***` is explicit:

> Beta-Carbon must be in an sp3 hybridization state for this feature to apply.

**Reference code.** It asks only that the beta-carbon carry no double bond. An aromatic
bond is not a `DOUBLE` bond, so a hydroxyl on an aromatic beta-carbon — a phenol —
satisfies the code's test and not the guidance's.

**What follows.** For `CN(N=O)c1ccccc1O` the guidance gives Category 2 and the code
gives Category 3. Phenolic NDSRIs are common, so this is not a corner case.

As with rule 3, **`balanced` follows the code despite the guidance being explicit**, and
for the same reason: agreement with the implementation that submissions are checked
against. `--rules guidance` applies the sp3 restriction as written. This reading is also
the more permissive one, since counting a phenol adds a deactivating score and raises
the limit.

---

## A fifth ambiguity that is not a divergence

The electron-withdrawing group list is the least precisely specified part of the CPCA,
but it is not a conflict between the two sources: Table B's footnote `**` excludes
carboxylic acid, aryl and ketone and then defers to Cross KP & Ponting DJ (2021),
*Comput Toxicol* **20**:100186, rather than enumerating the groups.

Because the guidance names no list, there is no guidance reading to switch to, so this
is not exposed as a setting. `cpca/ewg.py` reproduces the reference implementation's ten
patterns under every setting, including two quirks worth knowing about: FDA counts
matches rather than sides, so two groups on a single alpha-carbon score as "both sides";
and FDA computes an eleventh pattern, `ewg_aromaticBeta`, that it never adds to the
total. Both are reproduced so the two implementations agree.

Every EWG factor is marked `[needs review]` in the report for this reason. It is the one
Table B feature whose membership is a reading rather than a transcription.

## How the default is meant to be used

`balanced` is the default because it matches the implementation a submission is most
likely to be checked against, while refusing to score structures the published table
does not cover. It is not a claim that this is the correct reading of the guidance.

For strict conformance to the document as published, use `--rules guidance`. To
reproduce FDA's program exactly, including where it departs from the document, use
`--rules reference`. To see whether the choice changes a particular structure:

```
python score.py "<SMILES>" --compare-rules
```

Across the catalogue examined in this work, the settings differ on roughly one structure
in eighty.
