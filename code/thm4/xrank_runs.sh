#!/bin/sh
# driver for xrank.py (pilot_F2 step 1 falsifier); output runs/xrank.out.  Run from this folder.
O=runs/xrank.out
: > $O
for P in 2147483647 2147483629; do
 python xrank.py $P 1000 55 1 "5,7:R,R,R,R" >> $O
 python xrank.py $P 1000 55 2 "200001,200005:R,R,R,R" >> $O
 for S in "5,7:R,R,R,R;7,9:R,R,R,R" "5,7:R,R,R,R;9,11:R,R,R,R" "5,7:R,R,R,R;5,9:R,R,R,R" \
          "200001,200003:R,R,R,R;200001,200005:R,R,R,R" "200001,200005:R,R,R,R;300001,300007:R,R,R,R" \
          "200001,200003:R,R,R,R;200003,200005:R,R,R,R" "5,7:R,1,0,0;7,9:R,0,0,1" "5,7:R,0,1,0;5,9:R,1,0,0" \
          "5,7:0,0,0,1;7,9:0,1,0,0" "5,7:3,5,7,11;7,5:3,7,5,11" "5,7:3,5,7,11;-5,-7:3,-5,-7,11" "5,7:3,5,7,11;5,7:3,5,7,11" \
          "5,7:3,5,7,11;5,7:2,1,4,9"; do
  python xrank.py $P 1000 35 3 "$S" >> $O
 done
 for S in "5,7:R,R,R,R;7,9:R,R,R,R;5,9:R,R,R,R" "200001,200003:R,R,R,R;200003,200005:R,R,R,R;200001,200005:R,R,R,R" \
          "5,7:R,R,R,R;9,11:R,R,R,R;13,15:R,R,R,R" "5,7:3,5,7,11;7,9:R,R,R,R;7,5:3,7,5,11"; do
  python xrank.py $P 1000 25 4 "$S" >> $O
 done
done
echo DONE >> $O
