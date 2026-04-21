import sys
sys.path.insert(1, '/home/bzbauer/Documents/UArizona/HEP/Scripts/susypy')
import susypy as ss
import numpy as np
import matplotlib.pyplot as plt


input_name = "Bino-like_NLSP.txt"

susy_x = "/home/bzbauer/Programs/softsusy-4.1.22/softpoint.x"

input_dir = "/home/bzbauer/Documents/UArizona/HEP/Scripts/BinoNLSP/input"
output_dir = "/home/bzbauer/Documents/UArizona/HEP/Scripts/BinoNLSP/output"

slha = ss.SLHA(input_name,
                      susy_x,
                      in_dir=input_dir,
                      out_dir=output_dir)

masses = np.linspace(200, 450, 20)

mass_scan = ss.scan_params(slha, [("EXTPAR", 1)], [[f"{m:.8E}" for m in masses]])

particles = ["1000022", "1000023", "1000025", "1000035", "1000024", "1000037"]
labels = ["N1", "N2", "N3", "N4", "C1", "C2"]


# fig, ax = plt.subplots()

# ss.plot_scan(mass_scan, "EXTPAR", "1", "MASS", particles, abs_val=True, fig=fig, ax=ax, label_list=labels)

# ax.set_xlabel("Bino Mass [GeV]")
# ax.set_ylabel("Mass [GeV]")
# ax.set_title("Particle Masses vs. NLSP Mass")
# ax.legend()
# ax.grid()
# fig.savefig("electroweakinos_vs_bino_mass.png", dpi=300)

m1s = []
sigmas = [[],[],[]]
unctys = [[],[],[]]

for mass in mass_scan:
    m = np.abs(ss.gather_data([mass], "EXTPAR", "1"))
    sigma_n3n4, uncty_n3n4 = mass.cross_section("1000025", "1000035")
    sigma_c2c2, uncty_c2c2 = mass.cross_section("1000037", "-1000037")
    sigma_n3c2, uncty_n3c2 = mass.cross_section("1000025", "1000037")
    sigma_n4c2, uncty_n4c2 = mass.cross_section("1000035", "1000037")

    m1s.append(m)
    sigmas[0].append(sigma_n3n4 + sigma_n3c2)
    sigmas[1].append(sigma_c2c2 + sigma_n3c2 + sigma_n4c2)
    sigmas[2].append(sigma_n3n4 + sigma_n4c2)

    unctys[0].append(uncty_n3n4 + uncty_n3c2)
    unctys[1].append(uncty_c2c2 + uncty_n3c2 + uncty_n4c2)
    unctys[2].append(uncty_n3n4 + uncty_n4c2)


labels = ["N3", "C2", "N4"]

fig, ax = plt.subplots()
for sigma, uncty, label in zip(sigmas, unctys, labels):
    ax.errorbar(m1s, sigma, yerr=uncty, label=label)

ax.set_xlabel("Bino Mass [GeV]")
ax.set_ylabel("Cross Section [pb]")
ax.set_title("Particle Pair Production Cross Section vs Bino Mass")
ax.legend()
ax.grid()

plt.show()

