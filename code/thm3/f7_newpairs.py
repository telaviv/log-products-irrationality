"""RR-F7: exact ranges of new pairs (N, N+d), N >= 1, covered by Theorem C' (|m/n - 1| < 3e-6) but not by CDT's Theorem C
(< 1e-6), over all four representations (m, n) of the same product (negative pairs: log(1 + 1/(-k)) = -log(1 + 1/(k-1))).
Referee's newpairs.py with a wider search window (it checks that the window edges are not covered)."""
def cov(N, M, c):
    # |m/n - 1| < c * 1e-6  <=>  1e6 |m - n| < c |n|   (exact integers)
    reps = [(N, M), (M, N), (-(N + 1), -(M + 1)), (-(M + 1), -(N + 1))]
    return any(10 ** 6 * abs(m - n) < c * abs(n) for m, n in reps)
for d in (1, 2, 3, 10):
    lo = max(1, int(d / 3e-6) - 1000); hi = int(d / 1e-6) + 1000
    new = [N for N in range(lo, hi) if cov(N, N + d, 3) and not cov(N, N + d, 1)]
    assert new == list(range(new[0], new[-1] + 1)) and new[0] > lo and new[-1] < hi - 1
    print(f'd={d}: new N (pair N, N+{d}) = {new[0]:,} .. {new[-1]:,}  ({len(new)} values)')
