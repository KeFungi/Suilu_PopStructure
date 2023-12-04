from core_model import *
import sys
pop1=sys.argv[1]
pop2=sys.argv[2]

if len(sys.argv) >=5:
  omethod=sys.argv[3]
  imethod=sys.argv[4]

elif len(sys.argv) >=4:
  omethod=sys.argv[3]
  imethod=sys.argv[3]
  
else :
  omethod="tnc"
  imethod="tnc"

CEu_model.add_size_param("N0", N0=N0_neutral)
CEu_model.add_leaf("CEu-Belgium", N="N0", g=0)
CEu_model.optimize()

fit_model = CEu_model.copy()
fit_model.add_size_param("Ni1", N0=N0_intro)
fit_model.add_time_param("ti1", t0=t_intro)
fit_model.add_growth_param("g1", g0=0, lower=0, upper=0.1)
fit_model.add_leaf(pop1, N="Ni1", g="g1")
fit_model.move_lineages(pop1, "CEu-Belgium", t="ti1")

run_optimize(fit_model, method=imethod)

fit_model.add_size_param("Ni2", N0=N0_intro)
fit_model.add_time_param("ti2", t0=t_intro)
fit_model.add_growth_param("g2", g0=0, lower=0, upper=0.1)
fit_model.add_leaf(pop2, N="Ni2", g="g2")
fit_model.move_lineages(pop2, "CEu-Belgium", t="ti2")

run_optimize(fit_model, method=imethod)

fit_model.add_time_param("tp", t0=min(fit_model.get_params().ti1, fit_model.get_params().ti2)/2, upper_constraints=["ti1","ti2"])
fit_model.add_pulse_param("p12", p0=0.05, upper=0.5)
fit_model.move_lineages(pop2, pop1, t="tp", p="p12")

run_optimize(fit_model, method=omethod)

