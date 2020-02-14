#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phd.stats 
import phd.thermo
import phd.viz
import seaborn as sns
colors, palette = phd.viz.phd_style()
viridis = sns.color_palette('viridis', n_colors=4)

# Load the data and restrict
data = pd.read_csv('../../../data/RazoMejia2018_data.csv')
data = data[(data['repressors'] > 0) & (data['fold_change_A'] < 0.99) & 
        (data['fold_change_A']>0.01)]

# Compute the empirical bohr from each point 
data['empirical_F'] = -np.log(-1 + (1/data['fold_change_A'].values))
data['repressors'] *= 2

# Define the references
ref_r = 260
ref_c = 50 
ref_epR = -13.9

# Define the constants. 
ka = 139
ki = 0.53
ep_ai = 4.5


# Compute the reference and delta bohr for each 
bohr_c = phd.thermo.SimpleRepression(R=data['repressors'].values, 
                                    ep_r=data['binding_energy'].values,
                                    ka=ka, ki=ki, ep_ai=ep_ai,
                                    effector_conc=ref_c).bohr_parameter()
bohr_rep = phd.thermo.SimpleRepression(R=ref_r, 
                                    ep_r=data['binding_energy'].values,
                                    ka=ka, ki=ki, ep_ai=ep_ai,
                                    effector_conc=data['IPTG_uM'].values).bohr_parameter()
bohr_ep = phd.thermo.SimpleRepression(R=data['repressors'].values, 
                                    ep_r=ref_epR,
                                    ka=ka, ki=ki, ep_ai=ep_ai,
                                    effector_conc=data['IPTG_uM'].values).bohr_parameter()

# Compute the detlta
data['delF_c'] = data['empirical_F'] - bohr_c
data['delF_rep'] = data['empirical_F'] - bohr_rep
data['delF_ep'] = data['empirical_F'] - bohr_ep

# Group and compute the mean and sem. 
summarized = data.groupby(['repressors', 'operator', 
                           'IPTG_uM', 'binding_energy']).mean().reset_index()


#%%
# Set the colors for the strains
rep_colors = {22: colors['light_red'], 
              60: colors['light_brown'],  
              124: colors['light_green'], 
              260: colors['light_orange'], 
              1220: colors['light_purple'], 
              1740: colors['light_blue']} 

edge_colors = {22: colors['red'], 
              60: colors['brown'],  
              124: colors['green'], 
              260: colors['orange'], 
              1220: colors['purple'], 
              1740: colors['blue']} 

# Define the operators and their respective energies
glyphs = {'O1': 's', 'O2': '.', 'O3':'^'}

fig, ax = plt.subplots(1, 3, figsize=(7, 2.5), sharey=True)

c_plot = summarized.groupby(['operator', 'IPTG_uM']).agg(('mean', 'sem')).reset_index()
r_plot = summarized.groupby(['repressors', 'operator']).agg(('mean', 'sem')).reset_index()
ep_plot = summarized.groupby(['repressors', 'IPTG_uM']).agg(('mean', 'sem')).reset_index()

op_colors = {'O1':viridis[0], 'O2':viridis[1], 'O3':viridis[3]}
for g, d in c_plot.groupby(['operator']):
    ax[0].errorbar(np.log10(d['IPTG_uM'].values / ref_c), d['delF_c']['mean'], 
                   d['delF_c']['sem'], fmt='.', lw=1, color=op_colors[g], 
                   label=f'{g}', ms=6)

for g, d in r_plot.groupby(['operator']):
    ax[1].errorbar(np.log10(d['repressors'].values / ref_r), d['delF_rep']['mean'], 
                   d['delF_rep']['sem'], fmt='.', lw=1, color=op_colors[g], 
                   label=f'{g}', ms=6)
              
for g, d in r_plot.groupby(['repressors']):
    ax[2].errorbar(d['binding_energy']['mean'] - ref_epR, d['delF_ep']['mean'], 
                   d['delF_ep']['sem'], fmt='.', lw=1, color=edge_colors[g], 
                   label=f'{g}', ms=6)
            
    
# Plot the thoery curves
c_range = np.logspace(-2, 4, 200)
ep_range = np.linspace(-20, -8, 200)
rep_range = np.logspace(1, 3.5, 200)

ref = phd.thermo.MWC(ka=ka, ki=ki, ep_ai=ep_ai, effector_conc=ref_c).pact()
delta = phd.thermo.MWC(ka=ka, ki=ki, ep_ai=ep_ai, effector_conc=c_range).pact()
delF_c_theo = -np.log(delta / ref)
delF_rep_theo = -np.log(rep_range / ref_r)
delF_ep_theo = ep_range - ref_epR
ax[0].plot(np.log10(c_range / ref_c), delF_c_theo, 'k-', lw=1)
ax[1].plot(np.log10(rep_range / ref_r), delF_rep_theo, 'k-', lw=1)
ax[2].plot(ep_range - ref_epR, delF_ep_theo, 'k-', lw=1)

# Add legends. 
leg = ax[0].legend(loc='upper left', title='operator')
leg.get_title().set_fontsize(6)

leg = ax[1].legend(loc='upper right', title='operator')
leg.get_title().set_fontsize(6)

leg = ax[2].legend(loc='lower right', title='rep. per cell')
leg.get_title().set_fontsize(6)

# Add labels. 
ax[0].set_xlabel('$\log_{10} (c / c^{(ref)})$')
ax[0].set_ylabel('free energy shift [$k_BT$]')
ax[1].set_xlabel('$\log_{10} (R / R^{(ref)})$')
ax[2].set_xlabel(r'$\Delta\varepsilon_{RA} - \Delta\varepsilon_{RA}^{(ref)}$ [$k_BT$]')
plt.tight_layout()
plt.savefig('../figs/deltaF_reference.pdf', bbox_inches='tight')

# %%
# Set up the axes


# %%
