#!/bin/sh
# certified L13 continuation + rigorous determinant (F1's mono13.py / cert_det.py, unchanged copies), for the q given
cd "$(dirname "$0")"
for q in "$@"; do
  python -u mono13.py $q 2200 5 _cert 1700 diag > runs/mono13_q${q}_cert.log 2>&1
  python cert_det.py runs/mono13_q${q}_cert.pkl > runs/cert_det_q${q}.txt 2>&1
  echo "q=$q"; tail -n 1 runs/cert_det_q${q}.txt
done
