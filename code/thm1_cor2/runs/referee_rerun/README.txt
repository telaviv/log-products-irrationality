Gate referee's independent reruns of the L13 certificate (REFEREE_PILOT.md section 2.4), copied 2026-11-19 (RF1-F10)
from the referee scratch folder f1ref\ (scratch files may be cleared; these copies are the durable record).

Files
  mono13.py        referee-patched copy of the pilot's mono13_v1.py. Differences (env-controlled):
                     NOORTH=1 -> b = 2^-22 and tail bound (2n+1) 2^(20n) (4 log 2)/(q-1) WITHOUT the factor q^(-2n)
                                 (orthogonality-free initial data; N0 = 1248)
                     FR=<x>   -> Taylor step fraction (pilot 0.20)
                     NARC=<k> -> arcs per circle in circle_max (pilot 512)
                     plus asserts on every float conversion in circle_max (finite, nonzero).
  cert_det.py      identical to the pilot's cert_det.py.
  rec.py           identical to the pilot's rec.py (imported by mono13.py).
  floatrange.py    measures the float ranges of |A_R| and M_i on the circles (F-3).
  budget_ref.py    independent Arb recomputation of L, M (REFEREE_PILOT section 4).
  typeref.py       independent exact type check (REFEREE_PILOT section 3).
  ref_q32.log, mono13_q32_ref.pkl/.json, mono13_q32_ref_det.json, ref_det_q32.txt   -- q = 32 run
  ref_q26.log, mono13_q26_ref.pkl/.json, mono13_q26_ref_det.json, ref_det_q26.txt   -- q = 26 run

Settings of both runs: NOORTH=1 FR=0.16 NARC=768, i.e.
  NOORTH=1 FR=0.16 NARC=768 python -u mono13.py <q> 2200 5 _ref 1700 diag
  python cert_det.py runs/mono13_q<q>_ref.pkl
(the script reads runs/op_q<q>_25_13.json and runs/rec_q<q>.json relative to its own folder: to rerun, place it in the
pilot folder or copy those two files into a runs/ subfolder next to it).

Results (Arb, 2200 bits; the 14 greedy words '', 2, 3, A2, A3, A^2 2, A^2 3, A^3 2, A^3 3, A^4 2, A^4 3, A^5 2, A^6 2, A^7 2)
  q = 26: det = [257906186217.156323475542228360 +/- 3.7e-19] i   (same as the filed q = 26 certificate to 30 digits)
  q = 32: det = [992625.090850805422458198252898 +/- 3.83e-25] i, base 101, excludes 0
The in-script "CERTIFICATE" line in ref_q*.log is the default-word diagnostic (contains the dependent word A^5 3), not the
certificate; the certificate is *_det.json / ref_det_q*.txt (RF1-F4).
