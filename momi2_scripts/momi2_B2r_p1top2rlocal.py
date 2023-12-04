from core_model import *
import sys

pop1=sys.argv[1]
pop2=sys.argv[2]
logfile=sys.argv[3]

if len(sys.argv) >=6:
  omethod=sys.argv[4]
  imethod=sys.argv[5]

elif len(sys.argv) >=5:
  omethod=sys.argv[4]
  imethod=sys.argv[4]
  
else :
  omethod="tnc"
  imethod="tnc"

params=read_log(logfile)
CEu_model.add_size_param("N0", N0=params["N0"])
CEu_model.add_leaf("CEu-Belgium", N="N0", g=0)

fit_model = CEu_model.copy()
fit_model.add_size_param("Ni1", N0=params["Ni1"])
fit_model.add_time_param("ti1", t0=params["ti1"])
fit_model.add_growth_param("g1", g0=0, lower=0, upper=0.1)
fit_model.add_leaf(pop1, N="Ni1", g="g1")
fit_model.move_lineages(pop1, "CEu-Belgium", t="ti1")
fit_model.add_size_param("Ni2", N0=params["Ni2"])
fit_model.add_time_param("ti2", t0=params["ti2"], upper_constraints=["ti1"])
fit_model.add_growth_param("g2", g0=0, lower=0, upper=0.1)
fit_model.add_leaf(pop2, N="Ni2", g="g2")
fit_model.move_lineages(pop2, pop1, t="ti2")

run_optimize(fit_model, method=omethod)

