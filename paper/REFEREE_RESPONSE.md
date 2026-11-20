# Response to the final referee (REFEREE.md, accept after minor revision), 2026-11-19

## Part A: the stopped referee's edits

All twelve hunks are kept:
* A3 is amended: τ♯ also depends on 𝔪.
* A10a is corrected: the S_k are the numerators for the basis vectors; CDT's (14.3.6) has the single S₄₂.
* A4 and A9 are retypeset, to remove the four overfull boxes.

`check_y1.py` is kept, cited in §3.4 and listed in the README. `referee_partial_edits.patch` and `main_v1.tex` are deleted
(the template repositories keep no such copies).

## Part B: the findings

| # | point | action |
|---|---|---|
| F-1 | Most permissive representative | §1.1 now reads "(m, n) = (−N−1, −M−1), with \|m/n − 1\| = d/(M+1)". |
| F-2 | Bounds rounded the wrong way | Pooled table: 7.326965 → 7.326966, ≥ 0.01248 → ≥ 0.01247. §5.4: ≤ 16.95199 → ≤ 16.95200, numerator radius 1.5 → 1.6·10⁻¹⁵, CDT-map margin ≥ 0.031 → ≥ 0.0308. §4.1: M₁₄(26) ± 4 → ± 5·10⁻¹². Remark 4.6: ± 3.8 → ± 3.9·10⁻²⁵. README matched. |
| F-3 | Claims stated more widely than checked | §3.4: real, positive and distinct roots, and y₁ < 2⁻²⁰q², now stated for 2 ≤ q ≤ 60 (Arb, `check_y1.py`). A parenthesis says that for every q ≥ 2 only distinctness is proved, and that reality is used only for 26 ≤ q ≤ 53. §3.5 now fixes 2 ≤ q ≤ 60. Prop. 4.3 now speaks of "each of the four roots y_k of c₃₄". |
| F-4 | S_k wording | Fixed in §5.5 and in the proof of Prop. 5.2. |
| F-5 | "unpublished work" | Now "Apart from the preprints cited above (CDT, Lys18 and, for Lemma 6.2, LSS20), nothing in this paper depends on unpublished work by others." |
| F-6 | Reproducibility | §7 now says where each independent continuation rerun lives: q = 26, 32 in `code/thm1_cor2/runs/referee_rerun/`; q = 40, 53 in the gate referee's reruns, outputs not included. README changes: the `runcert.sh`/`mono13_v1.py` row corrected; rows added for `check_y1.py` and the discriminant; a note on the SymPy 1.14 failure in the last block of `indep_exp.py`; a note that certificate scripts overwrite filed outputs when re-run in place; the provenance paragraph now lists the one added file. |
| F-7 | Build warnings | `\relax` before the bracketed ball in Remark 4.6. |

## Build

`sh build.sh`:
* no errors;
* no undefined references;
* bibtex: no warnings;
* 0 overfull boxes;
* no LaTeX warnings.
