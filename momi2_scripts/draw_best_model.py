import momi
from matplotlib import pyplot as plt

base_N=1000
t_intro=200
N0_neutral=3000
N0_intro=500
mu_rate=1e-7
gen_time=1

def cleandraw(demoplot):
	demoplot.ax.xaxis.set_tick_params(labelsize=16)
	demoplot.ax.yaxis.set_tick_params(labelsize=12)
	demoplot.draw(rad=-0.2, alpha=1.0)
	demoplot.ax.set_ylabel('Time (years ago)', fontsize=16)
	demoplot.cbar.remove()
	demoplot.ax.set_xmargin(0.2)
	demoplot.ax.set_ymargin(0.1)
	#for lab in demoplot.ax.xaxis.get_majorticklabels():
	#	lab.set_ha("left")
	if len(demoplot.ax.texts) > 0:
		demoplot.ax.texts[0].set_size(24)
		demoplot.ax.texts[0].arrow_patch.set_lw(1.5)
		demoplot.ax.texts[1].set_size(12)
	plt.subplots_adjust(bottom=0.2, right=0.9, top=0.95, left=0.17)

#A1 SAm
NEu=33247.32568
Nin1=36649.82206
tin1=727.9621254
g1=0.004912346

A1 = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)

A1.add_leaf("Europe\n(Belgium)", N=NEu, g=0)
A1.add_leaf("South America\n(Bariloche)", N=Nin1, g=g1)
A1.move_lineages("South America\n(Bariloche)", "Europe\n(Belgium)", t=tin1)

A1plot =momi.DemographyPlot(
	A1,
    ["Europe\n(Belgium)", "South America\n(Bariloche)"],
    figsize=(6,8),
    draw=False,
    major_yticks=[100, 250, 500, 750, 1000]
    )
    
A1plot.base_N=base_N
cleandraw(A1plot)
plt.savefig("momi2_best_A1SAm", dpi=300)

#A1 NAm
NEu=33450.4866
Nin1=27864.23074
tin1=727.1308715
g1=0

A1 = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)

A1.add_leaf("Europe\n(Belgium)", N=NEu, g=0)
A1.add_leaf("North America\n(MA)", N=Nin1, g=g1)
A1.move_lineages("North America\n(MA)", "Europe\n(Belgium)", t=tin1)

A1plot =momi.DemographyPlot(
	A1,
    ["Europe\n(Belgium)", "North America\n(MA)"],
    figsize=(6,8),
    draw=False,
    major_yticks=[100, 250, 500, 750, 1000]
    )

A1plot.base_N=base_N
cleandraw(A1plot)
plt.savefig("momi2_best_A1NAm", dpi=300)


#B3 NZAU
NEu=32577.3621
Nin1=4562.453916
tin1=685.1519459
g1=0.0
Nin2=9759.025379
tin2=703.5331666
g2=0.003262121
tm=274.5980023
m=0.373998401

B3 = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)

B3.add_leaf("Europe\n(Belgium)", N=NEu, g=0)
B3.add_leaf("New Zealand", N=Nin1, g=g1)
B3.add_leaf("Australia\n(other)", N=Nin2, g=g2)
B3.move_lineages("New Zealand", "Europe\n(Belgium)", t=tin1)
B3.move_lineages("Australia\n(other)", "Europe\n(Belgium)", t=tin2)
B3.move_lineages("Australia\n(other)", "New Zealand", t=tm, p=m)


B3plot =momi.DemographyPlot(
	B3,
    ["Europe\n(Belgium)", "New Zealand", "Australia\n(other)"],
    figsize=(6,8),
    draw=False,
    major_yticks=[100, 250, 500, 750, 1000]
    )
B3plot.base_N=base_N
cleandraw(B3plot)
plt.savefig("momi2_best_B3AUNZ", dpi=300)



#B2 othertoWA
NEu=32931.27241
Nin1=17187.19345
tin1=747.6832452
gin1=0.002748821
Nin2=28596.31413
tin2=89.00728475
gin2=0.072879943

B2othertoWA = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)

B2othertoWA.add_leaf("Europe\n(Belgium)", N=NEu, g=0)
B2othertoWA.add_leaf("Australia\n(other)", N=Nin1, g=gin1)
B2othertoWA.add_leaf("Australia\n(WA)", N=Nin2, g=gin2)
B2othertoWA.move_lineages("Australia\n(other)", "Europe\n(Belgium)", t=tin1)
B2othertoWA.move_lineages("Australia\n(WA)", "Australia\n(other)", t=tin2)

B2othertoWAplot =momi.DemographyPlot(
	B2othertoWA,
    ["Europe\n(Belgium)", "Australia\n(other)", "Australia\n(WA)"],
    figsize=(6,8),
    draw=False,
    major_yticks=[100, 250, 500, 750, 1000]
    )
    
B3plot.base_N=base_N
cleandraw(B2othertoWAplot)
plt.savefig("momi2_best_B2othertoWA", dpi=300)


