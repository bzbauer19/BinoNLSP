import sys
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cplt
sys.path.insert(1, '/home/bzbauer/Documents/UArizona/HEP/Scripts/susypy')
import susypy as ss

dtime = dt.datetime.now()
date_string = (dtime.date().isoformat() 
    + "_" + str(dtime.time().hour) 
    + "-" + str(dtime.time().minute))

# Options

do_sigma = False

# Plotting palette definition

#color_palette = ["#648FFF", "#785EF0", "#DC267F",
#                "#FE6100", "#FFB000"]
                
#cmap_linear = cplt.LinearSegmentedColormap.from_list(
#    name="colorblind-ibm", colors=color_palette
#)

plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.Dark2.colors)

# Defining directories and file paths

input_name = "Bino-like_NLSP.txt"

susy_x = "/home/bzbauer/Programs/softsusy-4.1.22/softpoint.x"

input_dir = "/home/bzbauer/Documents/UArizona/HEP/Scripts/BinoNLSP/input"
output_dir = "/home/bzbauer/Documents/UArizona/HEP/Scripts/BinoNLSP/output"

plot_dir = "/home/bzbauer/Documents/UArizona/HEP/Scripts/BinoNLSP/Plots"

# Read in the base SLHA file

slha = ss.SLHA(input_name,
                      susy_x,
                      in_dir=input_dir,
                      out_dir=output_dir)

# Scan through a range of mu

masses = np.linspace(450, 800, 20)
tanbs = np.linspace(10, 60, 50)

mass_scan = ss.scan_params(slha, ("EXTPAR", 23), [f"{m:.8E}" for m in masses])
tanb_scan = ss.scan_params(slha, ("MINPAR", 3), [f"{tb:.8E}" for tb in tanbs])

# Define the electroweakino PDGs and labels

electroweakinos = ["1000022", "1000023", "1000025", "1000035", "1000024", "1000037"]
labels_ew = ["N1", "N2", "N3", "N4", "C1", "C2"]

# Plot electroweakino masses as a function of bino mass

fig, ax = plt.subplots()

fig, ax = ss.plot_scan(mass_scan, "EXTPAR", "23", "MASS", electroweakinos, abs_val=True, fig=fig, ax=ax, label_list=labels_ew)

ax.set_xlabel("Higgsino Mass [GeV]")
ax.set_ylabel("Mass [GeV]")
ax.set_title("Particle Masses vs. NLSP Mass")
ax.legend()
ax.grid()
fig.savefig(plot_dir + "/" + date_string + "electroweakinos_vs_higgsino_mass.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# Plot electroweakino masses as a function of bino mass

fig, ax = plt.subplots()

fig, ax = ss.plot_scan(tanb_scan, "MINPAR", "3", "MASS", electroweakinos, abs_val=True, fig=fig, ax=ax, label_list=labels_ew)

ax.set_xlabel(r"tan$\beta$")
ax.set_ylabel("Mass [GeV]")
ax.set_title("Particle Masses vs. NLSP Mass")
ax.legend()
ax.grid()
fig.savefig(plot_dir + "/" + date_string + "electroweakinos_vs_tanb.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# Compute the production cross sections for electroweakinos at LHC for bino NLSP
if do_sigma:
    m1s = []
    sigmas = [[],[],[],[],[],[]]
    unctys = [[],[],[],[],[],[]]

    for mass in mass_scan:
        m = np.abs(ss.gather_data([mass], "EXTPAR", "23"))
        sigma_n3n4, uncty_n3n4 = mass.cross_section("1000025", "1000035")
        sigma_c2c2, uncty_c2c2 = mass.cross_section("1000037", "-1000037")
        sigma_n3c2, uncty_n3c2 = mass.cross_section("1000025", "1000037")
        sigma_n4c2, uncty_n4c2 = mass.cross_section("1000035", "1000037")
        sigma_n3c2n, uncty_n3c2n = mass.cross_section("1000025", "-1000037")
        sigma_n4c2n, uncty_n4c2n = mass.cross_section("1000035", "-1000037")

        # Need help understanding the total cross section that we are looking for

        m1s.append(m)
        sigmas[0].append(sigma_n3n4)
        sigmas[1].append(sigma_c2c2)
        sigmas[2].append(sigma_n3c2)
        sigmas[3].append(sigma_n4c2)
        sigmas[4].append(sigma_n3c2n)
        sigmas[5].append(sigma_n4c2n)

        unctys[0].append(uncty_n3n4)
        unctys[1].append(uncty_c2c2)
        unctys[2].append(uncty_n3c2)
        unctys[3].append(uncty_n4c2)
        unctys[4].append(uncty_n3c2n)
        unctys[5].append(uncty_n4c2n)


    labels_sig = ["N3 + N4", "C2+ + C2-", "N3 + C2+", "N4 + C2+", "N3 + C2-", "N4 + C2-"]

    # Plot the total cross section for each  higgsino-like electroweakino 

    fig, ax = plt.subplots()
    #ax.set_prop_cycle('color', color_palette)
    for sigma, uncty, label in zip(sigmas, unctys, labels_sig):
        ax.errorbar(m1s, sigma, yerr=uncty, label=label)

    ax.set_xlabel("Higgsino Mass [GeV]")
    ax.set_ylabel("Cross Section [pb]")
    ax.set_title("Particle Pair Production Cross Section vs Higgsino Mass")
    ax.legend()
    ax.grid()

    fig.savefig(plot_dir + "/" + date_string + "electroweakino_sigma_vs_higgsino_mass.png", dpi=300, bbox_inches="tight")

# Branching Ratios maybe include the decay width?

for ew, name in zip(electroweakinos[1:], labels_ew[1:]):

    fig, ax = plt.subplots()

    fig, ax = ss.plot_scan(mass_scan, "EXTPAR", "23", ew, "all", abs_val=True, fig=fig, ax=ax, col_y=0)
    ax.set_xlabel("Higgsino Mass [GeV]")
    ax.set_ylabel("Branching Ratio")
    ax.set_title(name + " Branching Ratio vs Higgsino Mass")
    ax.grid()
    ax.legend()
    fig.savefig(plot_dir + "/" + date_string + "_" + name + "_BR_vs_higgsino_mass.png", dpi=300, bbox_inches="tight")
