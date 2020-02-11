#%%
import numpy as np
import pandas as pd 
import pickle
import phd.viz 
import holoviews as hv
import altair as alt
import bokeh.io 
import bokeh.plotting
hv.extension('bokeh')
_ = phd.viz.altair_theme()
colors, palette = phd.viz.bokeh_theme()
DATA_DIR = '../../../data'
# %%
# Load the sampler information
with open(f'{DATA_DIR}/SI_I_O2_R260.pkl', 'rb') as f:
    unpickler = pickle.Unpickler(f)
    gauss_flatchain = unpickler.load()
    ka = np.exp(-gauss_flatchain[:, 0][::100])
    ki = np.exp(-gauss_flatchain[:, 1][::100])

# Load the experimental data
data = pd.read_csv(f'{DATA_DIR}/RazoMejia2018_data.csv', comment='#')

# Make correction for tetramer to dimer
data = data[(data['repressors'] > 0) & (data['operator'] != 'Oid')].copy()
data['repressors'] *= 2

# Rename so things make some sense. 
data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'},
            inplace=True)

# Compute the aggregate properties. 
summary = data.groupby(['operator', 'repressors', 
                        'IPTGuM'])['fold_change'].agg(('mean', 'sem')).reset_index()
summary.rename(columns={'mean':'fc_mean', 'sem':'fc_sem'}, inplace=True)

summary['fc_min'] = summary['fc_mean'].values - summary['fc_sem'].values
summary['fc_max'] = summary['fc_mean'].values + summary['fc_sem'].values

# Isolate the fit strain for easy plotting
fit_strain = summary[(summary['operator']=='O2') & (summary['repressors']==260)]

# %%

# Define the inducer concentration range. 
c_range = np.logspace(-2, 4, 200)

# Set the colors for the strains
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

# Define the operators and their respective energies
operators = ['O1', 'O2', 'O3']
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7}

#%%
sampling_df = pd.DataFrame(np.array([ka, ki]).T, 
                          columns=['ka', 'ki'])

# Set up the dataframe for the fold-change. 
fc_df = pd.DataFrame()

# Iterate through each operator, repressor, and IPTG concentration to calculate
# the fold-change
for op, op_en in energies.items():
    for r, _ in rep_colors.items():
        for i, c in enumerate(c_range):
            arch = phd.thermo.SimpleRepression(r, op_en, ka=ka, ki=ki, 
                                               ep_ai=4.5, effector_conc=c)
            
            fc_min, fc_max = phd.stats.compute_hpd(arch.fold_change(), 0.95)
            fc_df = fc_df.append({'fc_min': fc_min,
                                  'fc_max': fc_max,
                                  'repressors': r,
                                  'IPTGuM': c,
                                  'operator': op,
                                  'binding_energy': op_en}, ignore_index=True)



#%%

# ############################################################################## 
# FIT STRAIN ONLY
# ##############################################################################

# Set up the plot base for the fit strain. . 
fit_base = alt.Chart(fit_strain[fit_strain['IPTGuM'] > 0])

points = fit_base.mark_point(size=20).encode(
                        x = 'IPTGuM:Q',
                        y = alt.Y('fc_mean:Q', axis=alt.Axis(title='fold-change in gene expression')),
                        fill=alt.value(colors['light_orange']),
                        stroke=alt.value(colors['dark_orange']),
                        strokeWidth=alt.value(0.5)
                    )

errs = fit_base.mark_errorbar().encode(
            x='IPTGuM',
            y=alt.Y('fc_min:Q', axis={'title':'fold-change in gene expression'}),
            y2='fc_max:Q',
            color=alt.Color('repressors:O', 
            scale=alt.Scale(range=(list(rep_colors.values()))), legend=None),
    )

fit_plot = errs + points

row = alt.hconcat()
for g, d in fc_df.groupby('operator'):
    plot_base = alt.Chart(d, width=150, height=150)
    error_band = plot_base.mark_area(opacity=0.45).encode(
        x=alt.X('IPTGuM:Q', 
                axis=alt.Axis(title='IPTG concentration [µM]', tickCount=4),
                scale=alt.Scale(type='log')),
        y=alt.Y('fc_min:Q',
                axis=alt.Axis(title='fold-change in gene expression'),
                scale=alt.Scale(domain=[-0.01, 1.15])),
        y2='fc_max:Q',
        fill=alt.Color('repressors:O', 
                        scale=alt.Scale(range=list(rep_colors.values())),
                        legend=alt.Legend(title='repressors per cell')),
        strokeWidth = alt.value(0.5),
        stroke = alt.Color('repressors:O',
                          scale=alt.Scale(range=list(edge_colors.values())))

       ).properties(
           title=f'operator {g}'
       )
    if g == 'O2':
        error_band += fit_plot
    row |= error_band

row
row.save('../assets/induction_fitstrain.svg')

#%%
# ##############################################################################
# ALL DATA 
# ##############################################################################

row = alt.hconcat()
for g, d in fc_df.groupby('operator'):
    plot_base = alt.Chart(d, width=150, height=150)
    error_band = plot_base.mark_area(opacity=0.45).encode(
        x=alt.X('IPTGuM:Q', 
                axis=alt.Axis(title='IPTG concentration [µM]', tickCount=4),
                scale=alt.Scale(type='log')),
        y=alt.Y('fc_min:Q',
                axis=alt.Axis(title='fold-change in gene expression'),
                scale=alt.Scale(domain=[-0.01, 1.15])),
        y2='fc_max:Q',
        fill=alt.Color('repressors:O', 
                        scale=alt.Scale(range=list(rep_colors.values())),
                        legend=alt.Legend(title='repressors per cell')),
        strokeWidth = alt.value(0.5),
        stroke = alt.Color('repressors:O',
                          scale=alt.Scale(range=list(edge_colors.values())))

       ).properties(
           title=f'operator {g}'
       )

    # Plot the data.  
    _data = summary[(summary['operator']==g) & (summary['IPTGuM'] > 0)]
    base = alt.Chart(_data)
    points = base.mark_point(size=20).encode(
        x = 'IPTGuM', 
        y= alt.Y('fc_mean:Q', axis={'title':'fold-change in gene expression'}),
        fill=alt.Color('repressors:O', 
                scale=alt.Scale(range=(list(rep_colors.values()))),
                legend=None),
        stroke = alt.Color('repressors:O',
                          scale=alt.Scale(range=list(edge_colors.values()))),
        strokeWidth=alt.value(0.5)
        )
    errs = base.mark_errorbar().encode(
            x='IPTGuM',
            y=alt.Y('fc_min:Q', axis={'title':'fold-change in gene expression'}),
            y2='fc_max:Q',
            color=alt.Color('repressors:O', 
                scale=alt.Scale(range=(list(rep_colors.values()))), legend=None),
    )
    row |= error_band + points + errs


row
row.save('../assets/induction_allstrain.svg')



# %%
# Bokeh implementation:
plots = []
bokeh.io.output_file('../assets/induction_fitstrain.html')
for g, d in fc_df.groupby('operator'):
    p = bokeh.plotting.figure(width=150, height=150, sizing_mode='scale_both',
                              x_axis_type='log', x_axis_label = 'IPTG [µM]',
                              y_axis_label='fold-change',
                              title=f'operator {g}',
                              x_range=[1E-2, 1E4],
                              toolbar_location=None)
    for r, _d in d.groupby('repressors'):
        p.varea(x='IPTGuM', y1='fc_min', y2='fc_max',
                        source=_d, fill_color=rep_colors[r], color=colors['black'],
                        legend_label=str(int(r)))
        p.line(x='IPTGuM', y='fc_min',
                        source=_d,  color=rep_colors[r],
                        line_width=0.75)
        p.line(x='IPTGuM', y='fc_max',
                        source=_d, color=rep_colors[r],
                        line_width=0.75)
        if (g == 'O2') & (r == 260):
            p.circle(x='IPTGuM', y='mean', fill_color='white', 
                    color=colors['dark_orange'], source=fit_strain, size=8)
    plots.append(p)
plots[0].legend.location = 'top_left'
plots[0].legend.title = 'repressors per cell'
plots[1].legend.visible = False
plots[2].legend.visible = False
plots[0].legend.label_text_font_size = "0.75em"
plots[0].legend.title_text_font_size = "0.75em"
for p in plots:
    p.axis.axis_label_text_font_size = "1em"
    p.title.text_font_size = "1em"
    p.axis.major_label_text_font_size = "0.75em"

lay = bokeh.layouts.gridplot([plots], sizing_mode='scale_both')
bokeh.io.save(lay)



# %%
