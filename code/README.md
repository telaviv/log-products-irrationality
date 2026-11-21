# Reproducibility code for "Products of two logarithms and arithmetic holonomy bounds"

This folder holds every script and certificate output used in the paper (`../paper/main.tex`). There is one directory
per result:

| directory | paper | content |
|---|---|---|
| `thm1_cor2/` | Section 3, Theorem 2 (25 <= N <= 31) cross-check (Remark 4.6) | Lysov's forms: recurrences, contact orders, minimal operators for q = 26..32, certified continuation, rank-14 certificates |
| `thm1_cor2/nonuniv/` | Section 4.6, Theorem 2 for 21 <= N <= 24 (new in v2) | q = 22..25: recurrences, minimal operators, continuation and rank-14 certificates; the degree-8 map, its preimage table and the certified bound for T(phi) |
| `thm1_cor2/fourterm/` | Section 4, Theorem 1 | local exponents, Casoratian, budget, the pipeline for q = 33..53, irreducibility certificates for q = 26..53 |
| `thm3/` | Section 5, Theorem 3 | Lemma S check, map certificates (cert2/cert3/cert3d), map optimisers, new-pair ranges |
| `thm3/f1/` | Section 5.5 | the 10 x 10 determinant of CDT Lemma 14.3.1, its factorisation and the factor Phi_6 |
| `thm4/` | Section 6, Theorem 4 | Lemma D, pooled tau, falsifier, pooled map certificates |

`thm1_cor2/` keeps the layout of the original working folder, because the `fourterm/` scripts read the q = 26..32
pickles and operators from `../runs/`. Run every script from inside its own folder.

## Requirements

Python 3.12, python-flint 0.9.0, sympy 1.14, mpmath 1.3, numpy 2.5, scipy 1.18 (scipy and numpy only for the map
optimisers, which are not part of any proof). `pip install python-flint sympy mpmath numpy scipy`. The `.sh` drivers
need a POSIX shell (Git Bash works on Windows).

## Provenance and changes

All files are byte copies of the working files of 19 November 2026, with three exceptions, where an absolute path to the
author's working folder was replaced by a relative one: `thm4/reopt.py` (path to `../thm3`), `thm4/tau_pool.py` (imports
`f1_tau.py`, copied into `thm4/` from the project's `cdt/` folder), and `thm1_cor2/fourterm/referee_4t/indep_exp.py`.
One file is new: `thm1_cor2/fourterm/check_y1.py`, added by the final referee for the root check of Section 3.4 (the
gate referee of Theorem 1 had done the same check in scratch).
Certificate outputs (`.json`, `.out`, `.log`) are unchanged, except that in version 2 the absolute paths of the
original working folders recorded in some of them (for example `"src"` fields and warning lines) were replaced by the
placeholder `<working folder>`; they record the run, not a path to use, and no script reads them. `thm3/lib_obs.py` imports `cdt4_budget.py` from another
project folder; it is used only by the map optimisers and by `xcheck_points.py` (not part of any proof) and is not
included.

Every certificate was re-run bit for bit by a referee on a separate copy (see Section 7 of the paper). The certificate
scripts write their outputs into `runs/` under the filed names, so a re-run in place overwrites the filed output; to compare,
run on a copy of the folder (the referees did this). No file here
exceeds 20 MB, so nothing is excluded from the repository.

## Map from the paper's numbers to scripts and outputs

### Section 3 (Lysov's forms), folder `thm1_cor2/`

| paper | script (command) | output |
|---|---|---|
| Relations for c_{n,l}, a_0, a_2, alpha (Prop. 3.4 proof) | `python i1_relations.py` | `runs/i1_relations.json` |
| Creative-telescoping certificates G_1, G_2, boundary cases | `python i1_certs.py` | `runs/i1_certs.json` |
| Order-4 operator p_0..p_4 (Prop. 3.4), comparison with `rec_q<q>.json` | `python i1_opalg.py` | `runs/i1_L4_symbolic.json` |
| Characteristic polynomial = 2^59 (8w+1)(9w-1)(25w-16) chi_w(-lambda); factors of p_0, p_4 | `python i1_extras.py` | `runs/i1_extras.txt` |
| Rodrigues identities, f_3 moments, contact orders n <= 24 (Lemma 3.1), Prop. 3.4 for n <= 2 | `python i1_contact.py` | `runs/i1_contact.json` |
| Recurrence at fixed q (found modularly; now the specialisation of the proved one) | `python rec.py q` | `runs/rec_q<q>.json` |
| Minimal operator L_13 (order 25, n-degree 13), exact check n <= 60 | `python findop.py q 25 13 240 60` | `runs/op_q<q>_25_13.json` |
| I1' (Prop. 3.6): zero right remainder, no nonnegative integer root of p_4 | `python i1_rec25.py q` | `runs/i1_rec25_q<q>.json` |
| Certified continuation, generator matrices g_A, g_2, g_3 (Section 4.5), q = 26..32 | the filed runs were made with `python -u mono13_v1.py q 2200 5 _cert 1700 diag`; `sh runcert.sh 26 27 ...` now calls the revision `mono13.py` with the same arguments (determinant-neutral, see the q = 26 `_rf3` row below) | `runs/mono13_q<q>_cert.{pkl,json,log}` |
| Rank-14 determinants (Remark 4.6), e.g. [257906186217.156323475542228360 +/- 3.7e-19] i at q = 26, [992625.090850805422458198252898 +/- 3.8e-25] i at q = 32 | `python cert_det.py runs/mono13_q<q>_cert.pkl` | `runs/mono13_q<q>_cert_det.json`, `runs/cert_det_q<q>.txt` (15-digit display) |
| Independent reruns at q = 26, 32 (orthogonality-free initial data, b = 2^-22, other discretisation) | see `runs/referee_rerun/README.txt` | `runs/referee_rerun/mono13_q{26,32}_ref*.{pkl,json}`, `ref_det_q{26,32}.txt` |
| Determinant-neutrality of the `mono13.py` revision (q = 26) | `python -u mono13.py 26 2200 5 _rf3 1700 diag` | `runs/mono13_q26_rf3*` |
| Budget M = (1 - 1/m) L - tau, q = 26..32 | `python budget_cert.py 26 27 ... 32` | printed |
| Denominator check (Lemma 3.3), exact, selected n | `runs/referee_rerun/typeref.py` | printed |

How to read the logs (float underflow of radii below 1e-308, diagnostic determinant lines): `runs/LOG_NOTES.txt`. The
q = 32 base-100 files `runs/*old_base100*` are an earlier run kept for the record; the certificate is the base-101 run.

### Section 4 (Theorem 1), folder `thm1_cor2/fourterm/`

| paper | script (command) | output |
|---|---|---|
| Exponent identity beta_k = 1/2, resultant cross-check (Prop. 4.3); Casoratian and p_0 roots for q = 2..60 (Lemma 4.4) | `python exponents.py` | `exponents.json`, `exponents.log` |
| Independent checks of the above (referee), including disc_y(c_34) = const w^3 (w-1)^2 (8w+1)^8 (9w-1)^8 (25w-16)^8 (Prop. 4.3) | `python referee_4t/indep_exp.py`, `python referee_4t/indep_cas.py` | printed (with SymPy 1.14 the last block of `indep_exp.py`, a second Casoratian computation, stops with a `PolynomialError` after everything else is printed; the Casoratian is checked by `exponents.py` and `indep_cas.py`) |
| Section 3.4: for q = 2..60 the roots of chi_{q^2}(-lambda) are real, positive and distinct, and y_1 < 2^-20 q^2 (used in Cor. 3.2 and Prop. 3.5) | `python check_y1.py` (added by the final referee) | printed: max y_1/(2^-20 q^2) = 0.353 at q = 2 |
| Table 1: M_14(q), m_min, y_2 (Arb) for q = 26..53 | `python budget_4t.py` | `budget_4t.json` |
| Pipeline q = 33..53: recurrence (specialised symbolic L4), L_13, I1', certified continuation, irreducibility | `sh run4t.sh 33 34 ... 53` (calls `mkrec.py`, `findop.py`, `i1_rec25.py`, `mono13.py`, `irred.py`) | `runs/rec_q<q>.json`, `runs/op_q<q>_25_13.json`, `runs/i1_rec25_q<q>.json`, `runs/mono13_q<q>_cert.{pkl,json,log}` |
| Table 1: det (a), det (c) (Prop. 4.5), q = 26..32 from `../runs/mono13_q<q>_cert.pkl` | `python irred.py ../runs/mono13_q<q>_cert.pkl runs/irred_q<q>.json` | `runs/irred_q<q>.json` (true radii), `.log` (display radii only) |
| Same, q = 33..53 | (in `run4t.sh`) | `runs/irred_q<q>.json` |
| Irreducibility from the independent reruns (Remark 4.6) | `python irred.py ../runs/referee_rerun/mono13_q<q>_ref.pkl runs/irred_q<q>_ref.json` | `runs/irred_q{26,32}_ref.json` |
| Table 1 (collected) | `python collect.py` | `summary_4t.json`, `runs/table_4t.md` |
| Denominator measurements quoted in the novelty discussion (not used in proofs) | `referee_4t/dens.py`, `referee_4t/dens2.py` | printed |

Times: `mono13.py` takes about 40 minutes per q (single core); everything else in this folder takes seconds to minutes.

### Section 4.6 (Theorem 2 for 21 <= N <= 24, i.e. q = 22..25; new in v2), folder `thm1_cor2/nonuniv/`

`mono13.py`, `cert_det.py`, `findop.py`, `i1_rec25.py` and `rec.py` are byte copies of the files in `thm1_cor2/`.
`make_rec.py` and `unit1_exponents.py` read `../runs/i1_L4_symbolic.json` (the proved symbolic order-4 operator).

| paper | script (command) | output |
|---|---|---|
| Recurrence at q = 22..25 as the specialisation of the proved symbolic L4 (q = 26 reproduces `../runs/rec_q26.json` exactly) | `python make_rec.py 22 23 24 25` | `runs/rec_q<q>.json` |
| L_13 and I1' at q = 22..25 (Prop. 4.8 proof) | `sh run_findop.sh 22 23 24 25` (calls `findop.py q 25 13 230 60`, then `i1_rec25.py`) | `runs/findop_q<q>.log`, `runs/op_q<q>_25_13.json`, `runs/i1_rec25_q<q>.json`, `runs/i1_rec25.log` |
| Certified continuation and the rank-14 determinant (Prop. 4.8, Table 2) | `sh runcert.sh 22 23 24 25` (calls `mono13.py q 2200 5 _cert 1700 diag`, then `cert_det.py`) | `runs/mono13_q<q>_cert.{pkl,json,log}`, `runs/mono13_q<q>_cert_det.json`, `runs/cert_det_q<q>.txt` |
| Preimage table (Lemma 4.9), L_phi, certified upper bound for T(phi) - L_phi, margin M_phi(q) (Table 2), also q = 21 (negative) and q = 26 (calibration) | `python cert_map.py q 512 1e-7` | `runs/cert_map_q<q>.json`, `runs/cert_map_q<q>.log` |
| Float cross-checks of T(phi) - L_phi (Jensen quadrature, direct double sampling); not part of the proof | `python proto_map.py q` | printed |
| Exponents 1/2 at y_1..y_4 for every q (a second check of Prop. 4.3) | `python unit1_exponents.py` | `runs/unit1_exponents.json` |

The `diag` determinant printed at the end of each `mono13` log uses a default word list and is a diagnostic only (it may
contain 0); the certificate is the `cert_det.py` output. Times: `mono13.py` about 50 minutes per q on one core (longer on
a loaded machine), `cert_map.py` a few seconds per q, everything else seconds to minutes. Keep at most two `mono13.py`
runs at a time on a machine with little free memory: a first attempt with four parallel runs ran out of memory at
q = 22, 23 (FLINT allocation error) and was repeated.

### Section 5 (Theorem 3), folder `thm3/`

| paper | script (command) | output |
|---|---|---|
| Lemma 5.1 (Lemma S), numerical spot check at boundary pairs | `python lemmaS_check.py` | printed |
| Section 5.4: certificate at 3e-6 (partition form; L = 4.79837381276..., numerator [9.23115527305872 +/- 1.52e-15], m <= 16.95200, margin >= 0.0015378) | `python cert3.py runs/p4_C_eps3.0e-06.json 3/1000000 1/2048,99/256,215/256,1 61/200,231/200,1109/500 0.002 16384` | `runs/cert3_p4_C_eps3.0e-06_eps3.0e-06_l3.json` (dyadic map coefficients in `c_dyadic`), `runs/cert3_C30.out` |
| Section 5.3/5.4: disc condition (c'), max log abs phi <= 13801311154345/2^40 = 12.55222 < log Y_1 | `python cert3d.py` (same arguments) | `runs/cert3d_p4_C_eps3.0e-06_eps3.0e-06_l3.json`, `runs/cert3d_C30.out` |
| Cross-checks at 2e-6, 2.5e-6 (margins >= 0.01138, >= 0.00131) and CDT's map at 1e-6 (>= 0.0308) | `python cert2.py runs/p4_C_eps2.0e-06.json 1/500000 ...` etc. | `runs/cert2_p4_C_eps2.0e-06_eps2.0e-06.json`, `runs/cert2_p4_C_eps2.5e-06_eps2.5e-06.json`, `runs/cert2_cdtfit_K64_eps1.0e-06.json`, `runs/cert2_C20.out`, `runs/cert2_C25.out` |
| New pairs (1.1): d = 2: 666,664 <= N <= 1,999,997; d = 3: 999,997 <= N <= 2,999,996; d = 10: 3,333,323 <= N <= 9,999,989 | `python f7_newpairs.py` | `runs/f7_newpairs.out` |
| Map optimisation (not part of the proof) | `p2*.py`-style drivers are not included; libraries `lib_map.py`, `lib_opt.py`, `lib_range.py`, `lib_band.py`, `lib_partition.py`, `lib_refine.py` | `runs/p4_*.json` (maps), `runs/p4_*.log` |
| Reproduction of CDT's own budget (sanity) | `python p0_reproduce.py` | `runs/p0_reproduce.json` |
| Limits (Section 8): margins at 4e-6 and 1e-5 | optimiser logs | `runs/p4_*.out`, `runs/p5_*.out` [numerical] |

The argument of `cert2.py`/`cert3.py` is the map JSON, epsilon as a fraction, then (cert3) the radii r_k and division
points theta_k, then delta_max and the number M0 of initial intervals. Runtime of the 3e-6 certificate: about 8 minutes.

### Section 5.5 (Phi_6; Proposition 5.2), folder `thm3/f1/`

The proof of Proposition 5.2 (the determinant identity) uses the three scripts of the first block; the second block is the
earlier interpolation route and the factor analysis, kept as cross-checks and for Phi_6's explicit form.

| paper | script | output |
|---|---|---|
| Prop. 5.2, operator: candidate L~ in Q[a,b,y][D] (y-degrees 12, 16, 17, 18, 19; total degree <= 12 in (a, b)); [y^k] L~(G_A) = 0 identically for k < 63; rank 63 of the ansatz at (7, 5), so CDT's L = L~; c_3/c_4 = CDT's printed log-derivative; agreement with `solve_L` at 6 fresh points | `python f1_Lsym.py` (27 s) | `runs/Lsym.txt` (L~), `runs/Lsym.out` |
| Prop. 5.2, a priori bound: S_k polynomial (deg_y 36+k, then 42; coefficients of total degree <= 28); o = (16,16,18,18,20,20,22,22,24,24), E = 200; det(a, b) = det M'(t, b)/t^200; totdeg det M' <= 990, deg_t <= 685, deg_b <= 690 | `python f1_detsym.py` (13 s) | `runs/detsym.pkl` (M' entries), `runs/detsym.out` |
| Prop. 5.2, identity: det M'(t, b_0) = kappa t^114 Pi(t + b_0, b_0) on all 692 lines b_0 = -346..345 (691 needed) | `python f1_identity.py` (about 35 min; resumable with `b_start b_stop`) | `runs/identity.json`, `runs/identity.out` ("lines checked: 692  all equal: True") |
| Operator L of CDT Prop. 14.1.11 at rational (a, b), point by point (library used by the above and by the earlier route) | `lib_f1.py`, `lib_L.py`, `lib_det.py`, `lib_fac.py` | (library) |
| Earlier route: det(a, b_0) (a - b_0)^86 for ten b_0, degree 373, 420 interpolation points + 6 checks (correct, but not a proof by itself: the a priori bound needs 686 points per line) | `python f1_det_univ.py` | `runs/det_b<b0>.json`, `.log` |
| Factors of the discriminant and resultants (no degree-6 factor); misprints in CDT (14.2.2), (14.2.3) | `python f1_factors.py` | `runs/factors.txt` |
| Phi_6 on a 7 x 7 grid; explicit Phi_6 | `python f1_phi6_grid.py`, `python f1_fit_phi6.py` | `runs/phi6.txt`, `runs/phi6_grid.out` |
| The determinant identity at 12 + 48 further points | `python f1_check_points.py` | `runs/check_points.json`, `runs/check_points.out` |
| Phi_6 = 16 mod 32 on all odd residue classes; diagonal bound; no zeros with abs(a), abs(b) <= 300 | `python f1_phi6_nonvanish.py` | `runs/phi6_nonvanish.out` |
| Transcriptions of CDT's R_10 and R_4 | | `r10.txt`, `r4.txt` |

### Section 6 (Theorem 4), folder `thm4/`

| paper | script (command) | output |
|---|---|---|
| Lemma 6.1 (Lemma D): c_3 = -2c_4 - 4 pi c_4, R_4(4) = (16a^2 - 40ab + 16b^2 + 9)^2, transcription of R_4 | `python lemmaD.py` | printed |
| Falsifier: exact rank mod two primes of the union of series (30 instances nullity 0; controls 10(D+1)) | `sh xrank_runs.sh` (calls `xrank.py`); self-test `python xrank_selftest.py` | `runs/xrank.out` |
| Pooled tau: K = 2: 1106977681/262702440, K = 3: 29601442477/7082066992 (xi = 169/100) | `python tau_pool.py` (uses `f1_tau.py`) | `runs/tau_pool.json` |
| Lemma S spot check at 1e-5, 3e-5, 1e-4 | `python lemmaS2.py` | printed |
| Certificate 5e-5, K = 2 (m = 27): L = 4.5159123241..., numerator [7.81991390318720 +/- 2.9e-16], m <= 25.885, margin >= 0.0124777 | `python cert_pool.py runs/reopt_e5e-5_m27.json 1/20000 1 - 27 1106977681/262702440 0.002 16384` | `runs/certpool_reopt_e5e-5_m27_eps5.0e-05_l0_m27.json`, `runs/certpool_5e-5_K2_l0.out` |
| Certificate 1e-4, K = 3 (m = 37): L = 4.4006730265..., numerator [7.32696512349209 +/- 1.44e-15], m <= 33.169, margin >= 0.0228725 | `python cert_pool.py ../thm3/runs/p4_D_eps1.0e-04.json 1/10000 1 - 37 29601442477/7082066992 0.002 16384` | `runs/certpool_p4_D_eps1.0e-04_eps1.0e-04_l0_m37.json`, `runs/certpool_1e-4_K3_l0.out` |
| Certificate 1e-5, K = 2 (margin >= 0.166084), map dilated by 1 - 2^-14 | `python dilate_map.py ../thm3/runs/p4_B_eps1.0e-05.json 0.99993896484375 runs/map_B_eps1e-05_r16383.json`, then `python cert_pool.py runs/map_B_eps1e-05_r16383.json 1/100000 1 - 27 1106977681/262702440 0.002 16384` | `runs/certpool_map_B_eps1e-05_r16383_eps1.0e-05_l0_m27.json`, `runs/certpool_1e-5_K2_l0.out` |
| Certificate 3e-5, K = 2 (margin >= 0.05585) | `python cert_pool.py ../thm3/runs/p4_D_eps3.0e-05.json 3/100000 1 - 27 1106977681/262702440 0.002 16384` | `runs/certpool_p4_D_eps3.0e-05_eps3.0e-05_l0_m27.json`, `runs/certpool_3e-5_K2_l0.out` |
| Smoke test: cert_pool with m = 17 and CDT's tau reproduces the Theorem 3 certificate bit for bit | `python cert_pool.py ../thm3/runs/p4_C_eps3.0e-06.json 3/1000000 1/2048,99/256,215/256,1 61/200,231/200,1109/500 17 1032659/242760 0.002 16384` | `runs/certpool_p4_C_eps3.0e-06_eps3.0e-06_l3_m17.json`, `runs/repro_cert3d_3e-6.out` |
| Re-optimised maps (not part of the proof); limits at 6e-5, 7e-5, 1e-4 (Section 8) | `python reopt.py eps 27 1106977681/262702440 src tag` | `runs/reopt_e{5,6,7}e-5_m27.{json,log}`, `runs/reopt_e1e-4_m27.*` |

Each certificate JSON stores the dyadic map coefficients it certified (`c_dyadic`), so a certificate can be re-checked
without the optimiser. Runtime: about 4 minutes per certificate.
