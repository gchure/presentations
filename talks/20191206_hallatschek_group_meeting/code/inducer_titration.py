#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.viz
import phd.thermo
colors, palette = phd.viz.phd_style()

# Define the number of repressors 
r = 200

# Define the inducer concentration range
c_range = np.logspace(-2, 4, 400)

# Define the DNA binding energies. 
ep_range = np.array([-16, -14, -12, -10])

# Mesh and compute the fold-change. 
c, ep = np.meshgrid(c_range, ep_range)
arch = phd.thermo.SimpleRepression(R=r, ep_r=ep, ep_ai=4.5, effector_conc=c,
                                ka=200, ki=1)
fc = arch.fold_change()


# %%
# Instantiate the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=100)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('inducer concentration relative to $K_A$ ($c / K_A$)', fontsize=10)
ax.set_ylabel('fold-change in gene expression', fontsize=10)
ax.tick_params(labelsize=9)


# Plot the foldchange curves. 
ep_colors = [colors['dark_red'], colors['dark_orange'], colors['orange'], 
            colors['light_orange']]

for i in range(4):
    ax.plot(c_range, fc[i, :], lw=2, color=ep_colors[i], label=int(ep_range[i]))

ax.legend(title=r'$\Delta\varepsilon_{RA}$ [$k_BT$]', fontsize=9)
plt.savefig('../figs/inducer_titration_theory.pdf', bbox_inches='tight',
            facecolor=None)
# %%
