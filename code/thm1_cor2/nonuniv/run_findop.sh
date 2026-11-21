#!/bin/sh
# exact L13 (order-25 recurrence, n-degree 13) at each q, then I1' by right division (i1_rec25.py)
cd "$(dirname "$0")"
for q in "$@"; do
  python -u findop.py $q 25 13 230 60 > runs/findop_q$q.log 2>&1
  tail -n 1 runs/findop_q$q.log
done
python -u i1_rec25.py "$@" > runs/i1_rec25.log 2>&1
tail -n 8 runs/i1_rec25.log
