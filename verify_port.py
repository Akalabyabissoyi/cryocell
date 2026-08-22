import sys, time; sys.path.insert(0,'.')
from cryocell.model import Params, simulate

def row(tag, **kw):
    S,R,E = simulate(Params(**kw))
    print(f"{tag:42s} imm {R['S_imm']*100:3.0f}%  24h {R['S_24']*100:3.0f}%  "
          f"| live {R['liveFrac']*100:3.0f} apop {R['apopFrac']*100:3.0f} necr {R['necrFrac']*100:3.0f}")
    return R

print("=== validation vs literature (must match the JS engine) ===")
row("hMSC suspension, no inhibitor")
row("+ Y-27632 5 uM   [Heng: 39.8->48.5]", rock_drug="y27632", rock_conc=5, rock_when="both")
row("+ Y-27632 10 uM  [Heng: 48.4]",       rock_drug="y27632", rock_conc=10, rock_when="both")
row("+ Y-27632 100 uM [Heng: 36.0, worse]",rock_drug="y27632", rock_conc=100, rock_when="both")
row("+ z-VAD-fmk      [best single]",      zvad=0.85)
row("4 C post-thaw hold [Heng 2006]",      T_recover=4)
row("shear-conditioned  [Bissoyi 2016]",   adhesion="sheared")
row("no CPA",                              cpa_key="none")
row("100 C/min",                           CR=100)
row("0.1 C/min",                           CR=0.1)
row("direct dilution",                     dilution="direct")
row("slow warming 5 C/min",                WR=5)
row("60 min at 37 C",                      T_add=37, hold_min=60)
row("trehalose only",                      cpa_key="tre")

print("\n=== cooling-rate sweep: inverted U ===")
best=(0,0)
for e in range(-4,13):
    r = 0.1*(10**(e/4))
    R = simulate(Params(CR=r))[1]
    if R['S_24']>best[1]: best=(r,R['S_24'])
    print(f"  {r:8.2f} C/min  S24 {R['S_24']*100:3.0f}%  ice {R['D_iif']*100:3.0f}  osm {R['D_osm']*100:3.0f}")
print(f"  optimum {best[0]:.2f} C/min at {best[1]*100:.0f}%")

print("\n=== timestep sensitivity (Cui's instability charge) ===")
for ds in [0.5, 1.0, 2.0, 4.0]:
    t0=time.time(); R = simulate(Params(dt_scale=ds))[1]; el=time.time()-t0
    print(f"  dt x{ds:<4} S24 {R['S_24']*100:6.2f}%   {el*1000:5.0f} ms")
