"""collect.py -- table of the per-q certificates for Theorem 4T: I1' (i1_rec25), findop exact check, irreducibility
(irred.py: both 13 x 13 Arb determinants exclude 0), budget M(14) > 0.  Writes summary_4t.json.
q = 26..32 use the pilot's pickles/ops (../runs) for continuation, and fourterm/runs/irred_q<q>.json; q = 33..53 are
fully produced in fourterm/runs.  Usage: python collect.py"""
import json, os, re
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, 'runs'); UP = os.path.join(H, '..', 'runs')
bud = json.load(open(os.path.join(H, 'budget_4t.json')))
rows = {}
for q in range(26, 54):
    r = dict(q=q, M14=bud[str(q)]['M14'], m_needed=bud[str(q)]['m_needed'])
    src = UP if q <= 32 else R
    try:
        j = json.load(open(os.path.join(src, 'i1_rec25_q%d.json' % q))); r['I1prime'] = j['remainder_zero'] and not j['p4_nonneg_integer_roots']
    except Exception as e: r['I1prime'] = None
    try:
        o = json.load(open(os.path.join(src, 'op_q%d_25_13.json' % q))); r['op_exact_check_failures'] = o['failures']; r['op_verified_to'] = o['verified_to']
    except Exception: r['op_exact_check_failures'] = None
    try:
        i = json.load(open(os.path.join(R, 'irred_q%d.json' % q)))
        r['irreducible'] = i['irreducible_certified']; r['det_i'] = i['det_i']; r['det_ii'] = i['det_ii']; r['pickle'] = i['pickle']
    except Exception: r['irreducible'] = None
    try:
        lg = open(os.path.join(R if q > 32 else UP, 'mono13_q%d_cert.log' % q)).read()
        m = re.search(r'check gA v = v: max \|diff\| upper (\S+)', lg); r['gA_v_check'] = m.group(1) if m else None
        m = re.search(r'true y_k (\[[^\]]*\]), base (\d+)', lg); r['yk_base'] = (m.group(1), m.group(2)) if m else None
    except Exception: pass
    r['ALL_OK'] = bool(r.get('I1prime') and r.get('op_exact_check_failures') == 0 and r.get('irreducible'))
    rows[q] = r
    print(q, 'I1\'=%s op_fail=%s irred=%s M14=%s m_needed=%s %s' % (r['I1prime'], r['op_exact_check_failures'], r['irreducible'], r['M14'], r['m_needed'], 'OK' if r['ALL_OK'] else '--'))
json.dump(rows, open(os.path.join(H, 'summary_4t.json'), 'w'), indent=1)
print('all OK:', [q for q, r in rows.items() if r['ALL_OK']])
print('not OK:', [q for q, r in rows.items() if not r['ALL_OK']])
