#!/usr/bin/env python3

"""Methotrexate example: how the output lost by each delivery arm depends on the DHFR and TS/SHMT levels.

Run from the repository root; prints a table (placeholder kinetics, see the project README).
"""

import sys,itertools,numpy as np
sys.path.insert(0,'devtools')
import sweep_common as sc, release_sweep as rs
data=sc.load_input('projects/methotrexate_dhfr/input.yml','/tmp/dbg/mtx')
b,_,m=sc.build_model(data)
def set_e(label,val):
    for sp in m.core_species:
        if sp.label.lower()==label.lower(): sp.concentration=sp.initial_concentration=val
base={'DHFR':1e-4,'TS':5e-4,'SHMT':5e-4}
print("DHFRx TS/SHMTx | bolus  liposome  depot  (% dTMP lost, 72h)")
for tsx in (0.2,1,3):
  for dx in (0.3,1,3,10,30):
    set_e('DHFR',base['DHFR']*dx); set_e('TS',base['TS']*tsx); set_e('SHMT',base['SHMT']*tsx)
    _,rows=rs.run_sweep(m,b,'methotrexate',['dTMP'],0.002,1.0,10000.0,259200.0)
    r={x['arm']:x['output_lost_pct'] for x in rows}
    print(f"{dx:5g} {tsx:5g} | {r['free bolus']:6.1f} {r['liposome']:7.1f} {r['slow depot']:7.1f}",flush=True)
