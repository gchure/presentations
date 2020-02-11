#%%
import numpy as np
import pandas as pd 
import pickle
import phd.viz 
import phd.thermo
import altair as alt

colors, palette = phd.viz.altair_theme()
DATA_DIR = '../../../data'

# Load the experimental data
data = pd.read_csv(f'{DATA_DIR}/RazoMejia2018_data.csv', comment='#')

# Make correction for tetramer to dimer
data = data[(data['repressors'] > 0) & (data['operator'] != 'Oid')].copy()
data['repressors'] *= 2

# Rename so things make some sense. 
data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'},
            inplace=True)

# Compute the aggregate properties. 
summary = data.groupby(['operator', 'repressors', 'binding_energy', 
                        'IPTGuM'])['fold_change'].agg(('mean', 'sem')).reset_index()
summary.rename(columns={'mean':'fc_mean', 'sem':'fc_sem'}, inplace=True)

summary['fc_min'] = summary['fc_mean'].values - summary['fc_sem'].values
summary['fc_max'] = summary['fc_mean'].values + summary['fc_sem'].values

# Compute the bohr parameter
summary['bohr'] = phd.thermo.SimpleRepression(summary['repressors'].values,
            summary['binding_energy'].values, ka=139, ki=0.53, ep_ai=4.5, 
            effector_conc=summary['IPTGuM'].values).bohr_parameter()

# %%
rep_colors = {22: colors['light_red'], 
              60: colors['light_brown'],  
              124: colors['light_green'], 
              260: colors['light_orange'], 
              1220: colors['light_purple'], 
              1740: colors['light_blue']} 

edge_colors = {22: colors['dark_red'], 
              60: colors['dark_brown'],  
              124: colors['dark_green'], 
              260: colors['dark_orange'], 
              1220: colors['dark_purple'], 
              1740: colors['dark_blue']}
base = alt.Chart(summary, width=300, height=150)

errors = base.mark_errorbar().encode(
        x = 'bohr:Q',
        y = alt.X('fc_min:Q', axis={'title':'fold-change'}),
        y2 = 'fc_max:Q',
        color = alt.Color('repressors:O', scale={'range':list(rep_colors.values())},
            legend=None),
)
points = base.mark_point().encode(
        x = alt.X('bohr:Q', axis={'title':'free energy [kT]'}),
        y = alt.Y('fc_mean:Q', axis={'title':'fold-change'}),
        shape = 'operator:N',
        fill = alt.Color('repressors:O', scale={'range':list(rep_colors.values())},
                    ),
        stroke = alt.Color('repressors:O', scale={'range':list(edge_colors.values())},
            legend=None),
        strokeWidth=alt.value(0.75)
)

# Compute and plot the master curve. 
bohr_range = np.linspace(-10, 10, 200)
fc = (1 + np.exp(-bohr_range))**-1
_df = pd.DataFrame({'bohr':bohr_range, 'fc':fc})
curve = alt.Chart(_df).mark_line().encode(
        x='bohr:Q',
        y=alt.Y('fc:Q', axis={'title':'fold-change'}),
        strokeWidth=alt.value(2)
        )
agg = points + errors + curve

agg.save('../assets/data_collapse.svg')

# %%



# %%
