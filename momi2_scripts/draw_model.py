import momi
from matplotlib import pyplot as plt

t_intro=200
N0_neutral=3000
N0_intro=500
mu_rate=1e-6
gen_time=10

yticks = [100, 150, 200, 250, 300]
fsize=(6,8)

def cleandraw(demoplot):
	demoplot.ax.xaxis.set_tick_params(labelsize=20)
	demoplot.ax.yaxis.set_tick_params(labelsize=16)
	demoplot.draw()
	demoplot.ax.set_ylabel('Time (years ago)', fontsize=20)
	demoplot.cbar.remove()
	if len(demoplot.ax.texts) > 0:
		demoplot.ax.texts[0].set_size(28)
		demoplot.ax.texts[0].arrow_patch.set_lw(2)
		demoplot.ax.texts[1].set_size(16)
	plt.subplots_adjust(bottom=0.15, right=0.9, top=0.95, left=0.17)

A1 = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)
A1.add_leaf("Europe", N=N0_neutral, g=0)
A1.add_leaf("introduced", N=N0_intro, g=0.03)
A1.move_lineages("introduced", "Europe", t=t_intro)


A1plot =momi.DemographyPlot(
    A1,
    ["Europe", "introduced"],
    figsize=fsize,
    draw=False,
    major_yticks=yticks
    )
    
cleandraw(A1plot)
plt.savefig("momi2_A1", dpi=300)

B1 = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)
B1.add_leaf("Europe", N=N0_neutral, g=0)
B1.add_leaf("introduced1", N=N0_intro, g=0.03)
B1.add_leaf("introduced2", N=N0_intro, g=0.03)
B1.move_lineages("introduced1", "Europe", t=t_intro*0.8)
B1.move_lineages("introduced2", "Europe", t=t_intro)

B1plot =momi.DemographyPlot(
	B1,
    ["Europe", "introduced1", "introduced2"],
    figsize=fsize,
    draw=False,
    major_yticks=yticks
    )

cleandraw(B1plot)
plt.savefig("momi2_B1", dpi=300)
    
B2 = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)
B2.add_leaf("Europe", N=N0_neutral, g=0)
B2.add_leaf("introduced1", N=N0_intro, g=0.02)
B2.add_leaf("introduced2", N=N0_intro, g=0.02)
B2.move_lineages("introduced1", "Europe", t=t_intro)
B2.move_lineages("introduced2", "introduced1", t=t_intro*0.8)

B2plot =momi.DemographyPlot(
	B2,
    ["Europe", "introduced1", "introduced2"],
    figsize=fsize,
    draw=False,
    major_yticks=yticks
    )
    
cleandraw(B2plot)
plt.savefig("momi2_B2", dpi=300)

B3 = momi.DemographicModel(N_e=N0_neutral, gen_time=gen_time, muts_per_gen=mu_rate)
B3.add_leaf("Europe", N=N0_neutral, g=0)
B3.add_leaf("introduced1", N=N0_intro, g=0.03)
B3.add_leaf("introduced2", N=N0_intro, g=0.03)
B3.move_lineages("introduced1", "Europe", t=t_intro*0.8)
B3.move_lineages("introduced2", "Europe", t=t_intro)
B3.move_lineages("introduced2", "introduced1", t=t_intro*0.3, p=0.05)


B3plot =momi.DemographyPlot(
	B3,
    ["Europe", "introduced1", "introduced2"],
    figsize=fsize,
    draw=False,
    major_yticks=yticks
    )
    
cleandraw(B3plot)
plt.savefig("momi2_B3", dpi=300)