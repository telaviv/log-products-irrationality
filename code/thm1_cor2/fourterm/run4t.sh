#!/bin/sh
# per q: order-4 recurrence (specialised symbolic L4) -> minimal operator L13 (findop) -> I1' (i1_rec25) ->
# certified continuation (mono13, as runcert.sh) -> irreducibility certificate (irred.py)
cd "$(dirname "$0")"
for q in "$@"; do
  python mkrec.py $q > runs/mkrec_q$q.log 2>&1
  python findop.py $q 25 13 240 60 > runs/findop_q$q.log 2>&1
  python i1_rec25.py $q > runs/i1_rec25_q$q.log 2>&1
  python -u mono13.py $q 2200 5 _cert 1700 diag > runs/mono13_q${q}_cert.log 2>&1
  python irred.py runs/mono13_q${q}_cert.pkl runs/irred_q$q.json > runs/irred_q$q.log 2>&1
  echo "q=$q"; tail -n 1 runs/irred_q$q.log
done
