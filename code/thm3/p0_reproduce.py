"""p0_reproduce.py -- step 1: reproduce CDT4_BUDGET d=1 numbers at t = 1-1e-4 and 1-1e-6 with the existing scripts
(imported read-only from schanuel/five_exp/cdt4). Output runs/p0_reproduce.json. [Numerical]"""
import sys, os, json, time
B = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'schanuel', 'five_exp', 'cdt4'))
sys.path.insert(0, B)
import cdt4_budget as CB
out = {}
for lab, aC, bC in [('1-1e-4', 2 * 10**4 - 1, 2 * (10**4 - 1) - 1), ('1-1e-6', 2 * 10**6 - 1, 2 * (10**6 - 1) - 1)]:
    t0 = time.time()
    rec = CB.instance(lab, aC, bC, ds=(1,))
    b = rec['oc_explicit']; r1 = rec['budget'][1]
    print(lab, 'delta_oc', rec['delta_oc'], 'L_oc_hi', rec['L_oc_hi'], 'best', b['kind'], b['L'], b['T'],
          'M_expl', r1['M_explicit'], 'M_ceil', r1['M_ceiling'], 'm_max', r1['m_max'], f'{time.time()-t0:.0f}s', flush=True)
    out[lab] = dict(delta_oc=rec['delta_oc'], L_oc_hi=rec['L_oc_hi'], best=b, tried=rec['oc_tried'], budget1=r1, sing=rec['sing'])
json.dump(out, open('runs/p0_reproduce.json', 'w'), indent=1, default=str)
