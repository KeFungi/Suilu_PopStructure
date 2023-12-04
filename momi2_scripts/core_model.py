import momi

t_intro=200
N0_neutral=3000
N0_intro=500
mu_rate=1e-6
gen_time=10

sfs = momi.Sfs.load("sfs.gz")

CEu_model = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)
CEu_model.set_data(sfs)

def run_optimize(mymodel, **kwargs):
  opt_log = mymodel.optimize(options={'maxiter':1000}, **kwargs)
  print(opt_log, flush=True)
  print({'mu':mu_rate, 'gen':gen_time}, flush=True)
  print(mymodel.log_likelihood(), flush=True)
  print(mymodel.get_params(), flush=True)

def read_log(logfile):
	with open(logfile) as f:
	    for line in f:
	        pass
	    last_line = line
	    params=eval(last_line.rsplit("(")[1].rsplit(")")[0])
	    print(params, flush=True)
	    return(params)