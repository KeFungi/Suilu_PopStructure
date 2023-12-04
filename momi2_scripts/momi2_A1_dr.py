from core_model import *
import sys
pop=sys.argv[1]

if len(sys.argv) >=3:
  omethod=sys.argv[2]
else:
  omethod="tnc"

CEu_model.add_size_param("N0", N0=N0_neutral)
CEu_model.add_leaf("CEu-Belgium", N="N0", g=0)
CEu_model.optimize()

fit_model = CEu_model.copy()
fit_model.add_size_param("Ni", N0=N0_intro)
fit_model.add_time_param("ti", t0=t_intro)
fit_model.add_growth_param("g_in", g0=0, lower=0, upper=0.1)
fit_model.add_leaf(pop, N="Ni", g="g_in")
fit_model.move_lineages(pop, "CEu-Belgium", t="ti", g=0)

run_optimize(fit_model, method=omethod)
