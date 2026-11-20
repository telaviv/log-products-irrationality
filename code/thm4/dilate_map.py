"""dilate_map.py src.json r out.json -- the map phi_r(z) = phi(r z), i.e. psi_r(z) = z exp(log r + c_0 + sum_k c_k r^k z^k).
Used to move the transversal crossing of phi(T) with the positive real axis at y = 100000.06 (range-pilot map p4_B, 10^-5)
below Y1 = 1/eps - 6 = 99994.  Floats in, floats out; cert_pool.py rounds the output to dyadics and certifies that object."""
import sys, json, math
src, r, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
d = json.load(open(src))
c = [ck * r ** k for k, ck in enumerate(d['c'])]
c[0] += math.log(r)
json.dump(dict(eps=d['eps'], K=d['K'], c=c, src=src, dilation=r, L=math.log(256) + c[0]), open(out, 'w'), indent=1)
print('L =', math.log(256) + c[0])
