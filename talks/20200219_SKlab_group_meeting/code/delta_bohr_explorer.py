#%%
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import phd.viz
import phd.thermo
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, IndexFilter
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button
from bokeh.embed import components
phd.viz.bokeh_theme()
colors, palette = phd.viz.phd_style()
bokeh.plotting.output_file("../figs/delta_bohr_explorer.html", mode='inline')

# Define the invariant properties
ka = 200
ki = 1
ep_ai = 5
n_points = 200
c_range = np.logspace(-5, 5, n_points)
c_ki_range = c_range / ki

# Define the perturbed state
rep_slider = Slider(title='repressors per cell', start=1, end=1000, step=10, 
                    value=200, bar_color=colors['light_red'])
ep_slider = Slider(title='DNA binding energy [kT]', start=-25, end=0, 
                    step=0.5, value=-14, bar_color=colors['light_orange'])
c_slider = Slider(title='log\u2081\u2080 (c / Ki)', start=-4, end=6, step=0.01, value=0, bar_color=colors['light_blue'])


# Define the reference state
ref_rep_slider = Slider(title='reference repressors per cell', start=1, end=1000, step=10, 
                    value=200, bar_color=colors['black'])
ref_ep_slider = Slider(title='reference DNA binding energy [kT]', start=-25, end=0, 
                    step=0.5, value=-14, bar_color=colors['black'])
ref_c_slider = Slider(title='log\u2081\u2080 (c_ref / Ki)', 
        start=-4, end=6, step=0.01, value=0, bar_color=colors['black'])

# Define the architectures
ref_arch_point = phd.thermo.SimpleRepression(R=ref_rep_slider.value, ep_r=ref_ep_slider.value,
                                       ka=ka, ki=ki, ep_ai=ep_ai, 
                                       effector_conc=10**(ref_c_slider.value * ki))
ref_MWC = phd.thermo.MWC(ka=ka, ki=ki, ep_ai=ep_ai, effector_conc=ref_c_slider.value)

arch_point = phd.thermo.SimpleRepression(R=rep_slider.value, ep_r=ep_slider.value,
                                       ka=ka, ki=ki, ep_ai=ep_ai, 
                                       effector_conc=10**(c_slider.value * ki))
MWC = phd.thermo.MWC(ka=ka, ki=ki, ep_ai=ep_ai, effector_conc = c_range)
ref_arch_curve = phd.thermo.SimpleRepression(R=ref_rep_slider.value, ep_r=ref_ep_slider.value,
                                       ka=ka, ki=ki, ep_ai=ep_ai, 
                                       effector_conc=c_ki_range)
arch_curve = phd.thermo.SimpleRepression(R=rep_slider.value, ep_r=ep_slider.value,
                                       ka=ka, ki=ki, ep_ai=ep_ai, 
                                       effector_conc=c_ki_range)

# Define the source for the curves
ref_fc = ref_arch_curve.fold_change()
fc = arch_curve.fold_change()
induction_source = ColumnDataSource({'c':c_range, 
                                     'c_ki':np.log10(c_ki_range), 
                                     'fc':fc, 
                                     'fc_ref':ref_fc})

# Define the point sources.
ref_fc_point = ref_arch_point.fold_change()
fc_point = arch_point.fold_change()
ref_bohr = ref_arch_point.bohr_parameter()
bohr = arch_point.bohr_parameter()
delF = bohr - ref_bohr
point_source = ColumnDataSource({'c':[10**(c_slider.value * ki)],
                                 'c_ki':[c_slider.value],
                                 'c_ref':[10**(ref_c_slider.value * ki)],
                                 'c_ki_ref':[ref_c_slider.value],
                                 'c_c_ref':[np.log10((10**(c_slider.value * ki))/(10**(ref_c_slider.value * ki)))],
                                 'r':[rep_slider.value],
                                 'r_r_ref':[np.log10(rep_slider.value / ref_rep_slider.value)],
                                 'ep':[ep_slider.value],
                                 'ep_ep_ref':[ep_slider.value - ref_ep_slider.value],
                                 'fc':[fc_point],
                                 'fc_ref':[ref_fc_point],
                                 'bohr': [bohr],
                                 'bohr_ref':[ref_bohr],
                                 'delF_c': [delF],
                                 'delF_r': [delF],
                                 'delF_ep': [delF]})

# Define the delF curve sources. 
rep_range = np.logspace(0, 6,  n_points)
r_r0_range = np.log10(rep_range / ref_rep_slider.value)
c_c0_range = np.log10(c_range / 10**(ref_c_slider.value * ki))
ep_range = np.linspace(-25, 0,  n_points)
delF_ep = ep_range - ref_ep_slider.value 

delF_source = ColumnDataSource({'rep_range':rep_range, 
                                 'r_r_ref': np.log10(rep_range / ref_rep_slider.value),
                                 'c_range':c_range, 
                                 'c_c_ref':np.log10(c_range / 10**(ref_c_slider.value * ki)),
                                 'ep_range':ep_range,
                                 'ep_ep_ref': ep_range - ref_ep_slider.value,
                                 'delF_c':-np.log(MWC.pact() / ref_MWC.pact()),
                                 'delF_r':-np.log(rep_range / ref_rep_slider.value),
                                 'delF_ep':ep_range - ref_ep_slider.value})

# Set up the figure canvas
plot_width = 375
plot_height = 325
p_fc = bokeh.plotting.figure(width=plot_width, height=plot_height, 
                             x_axis_label='log\u2081\u2080 (c / Ki)',
                             y_axis_label = 'fold-change',
                             y_range=[-0.01, 1.1],
                             toolbar_location=None)

p_bohr = bokeh.plotting.figure(width=plot_width, height=plot_height, 
                             x_axis_label='free energy [kT]',
                             y_axis_label = 'fold-change',
                             toolbar_location=None)

p_dF_rep = bokeh.plotting.figure(width=plot_width, height=plot_height, 
                             x_axis_label='log\u2081\u2080 (R/R_ref)',
                             y_axis_label = 'free energy shift [kT]',
                             toolbar_location=None)

p_dF_c = bokeh.plotting.figure(width=plot_width, height=plot_height, 
                             x_axis_label='log\u2081\u2080(c/ c_ref)',
                             y_axis_label = 'free energy shift [kT]',
                             toolbar_location=None)

p_dF_ep = bokeh.plotting.figure(width=plot_width, height=plot_height, 
                             x_axis_label='∆ε - Δε_ref [kT]',
                             y_axis_label='free energy shift [kT]',
                             toolbar_location=None)

# Invariant master curve
bohr_range = np.linspace(-15, 15,  n_points)
master_curve = (1 + np.exp(-bohr_range))**-1
p_bohr.line(bohr_range, master_curve, color=colors['black'], line_width=2)

# delF curves 
p_dF_rep.line('r_r_ref', 'delF_r', color=colors['black'], line_width=2, source=delF_source)
p_dF_c.line('c_c_ref', 'delF_c', color=colors['black'], line_width=2, source=delF_source)
p_dF_ep.line('ep_ep_ref', 'delF_ep', color=colors['black'], line_width=2, source=delF_source)


# Foldchange curves. 
p_fc.line('c_ki', 'fc_ref', color=colors['black'], line_width=2, 
          source=induction_source, legend_label='reference')
p_fc.line('c_ki', 'fc', color=colors['purple'], line_width=2, 
          source=induction_source)
p_fc.circle('c_ki_ref', 'fc_ref', color='white', line_color=colors['black'],
            size=10, line_width=2, source=point_source)
p_fc.circle('c_ki', 'fc', color='white', line_color=colors['purple'],
            size=10, line_width=2, source=point_source)
p_fc.legend.location = 'top_left'

# Bohr points. 
p_bohr.circle('bohr_ref', 'fc_ref', color='white', line_color=colors['black'],
            size=10, line_width=2, source=point_source)
p_bohr.circle('bohr', 'fc', color='white', line_color=colors['purple'],
            size=10, line_width=2, source=point_source)

#delF points
p_dF_c.circle('c_c_ref', 'delF_c', color='white', line_color=colors['purple'],
            size=10, line_width=2, source=point_source)
p_dF_rep.circle('r_r_ref', 'delF_r', color='white', line_color=colors['purple'],
            size=10, line_width=2, source=point_source)
p_dF_ep.circle('ep_ep_ref', 'delF_ep', color='white', line_color=colors['purple'],
            size=10, line_width=2, source=point_source)


# Define the callback. 
js = """

// Define the parameter values
var r_ref = ref_rep_slider.value;
var r = rep_slider.value;
var c_ref = Math.pow(10, ref_c_slider.value * ki);
var c = Math.pow(10, c_slider.value * ki);
var ep_ref = ref_ep_slider.value;
var ep = ep_slider.value;


// Define functions for easy calculation of properties
function computePact(c, ka, ki, ep_ai, n) {
    var numer = Math.pow(1 + c / ki, n);
    var denom = Math.pow(1 + c / ka, n);
    return Math.pow(1 + Math.exp(-ep_ai) * numer / denom, -1);
    }

function computeFoldChange(R, ep_r, c, ka, ki, ep_ai, n) {
    var pact = computePact(c, ka, ki, ep_ai, n);
    // Note that 4.6E6 is hard coded as the number of nonspecific 
    // sites
    return Math.pow(1 + pact * (R / 4600000) * Math.exp(-ep_r), -1);
    }

function computeBohrParameter(R, ep_r, c, ka, ki, ep_ai, n) {
    var pact = computePact(c, ka, ki, ep_ai, n);
    // Note that 4.6E6 is hard coded as the number of nonspecific 
    // sites
    return -Math.log(pact) - Math.log(R/4600000) + ep_r;
    }

// Update the fold-change curves. 
for (var i = 0 ; i < ind_source.data['c'].length; i++) {
    ind_source.data['fc'][i] = computeFoldChange(r, ep, c_range[i], ka, ki, ep_ai, n);
    ind_source.data['fc_ref'][i] = computeFoldChange(r_ref, ep_ref, c_range[i], ka, ki, ep_ai, n);
}

// Update the collapse curves. 
for (var i = 0; i < delF_source.data['c_range'].length; i++ ) {
    delF_source.data['delF_c'][i] = -Math.log(computePact(delF_source.data['c_range'][i], ka, ki, ep_ai, n) / computePact(c_ref, ka, ki, ep_ai, n));
    delF_source.data['c_c_ref'][i] = Math.log10(delF_source.data['c_range'][i] / c_ref);
    delF_source.data['delF_r'][i] = -Math.log(delF_source.data['rep_range'][i] / r_ref);
    delF_source.data['r_r_ref'][i] = Math.log10(delF_source.data['rep_range'][i] / r_ref)
    delF_source.data['delF_ep'][i] = delF_source.data['ep_range'][i] - ep_ref;
    delF_source.data['ep_ep_ref'][i] = delF_source.data['ep_range'][i] - ep_ref;
    }
// Update the fold-change points
point_source.data['fc'] = [computeFoldChange(r, ep, c, ka, ki, ep_ai, n)];
point_source.data['c_ki'] = [Math.log10(c / ki)];
point_source.data['fc_ref'] = [computeFoldChange(r_ref, ep_ref, c_ref, ka, ki, ep_ai, n)];
point_source.data['c_ki_ref'] = [Math.log10(c_ref / ki)];

// Update the bohr points
point_source.data['bohr'] = [computeBohrParameter(r, ep, c, ka, ki, ep_ai, n)];
point_source.data['bohr_ref'] = [computeBohrParameter(r_ref, ep_ref, c_ref, ka, ki, ep_ai, n)];

// Compute the free energy shifts. 
point_source.data['delF_c'] = [-Math.log(computePact(c, ka, ki, ep_ai, n) / computePact(c_ref, ka, ki, ep_ai, n))];
point_source.data['c_c_ref'] = [Math.log10(c / c_ref)];
point_source.data['delF_r'] = [-Math.log(r / r_ref)];
point_source.data['r_r_ref'] = [Math.log10(r / r_ref)];
point_source.data['delF_ep'] = [ep - ep_ref];
point_source.data['ep_ep_ref'] = [ep - ep_ref];

// Emit the changes. 
ind_source.change.emit();
point_source.change.emit();
delF_source.change.emit();
"""

# Define arguments for the callback
args = {'ka':ka, 'ki':ki, 'ep_ai':ep_ai, 'n':2,
        'c_range':c_range, 'ind_source':induction_source,
        'point_source':point_source, 'delF_source':delF_source,
        'rep_slider':rep_slider, 'c_slider':c_slider, 'ep_slider':ep_slider,
        'ref_rep_slider':ref_rep_slider, 'ref_c_slider':ref_c_slider, 'ref_ep_slider':ref_ep_slider}

# Define the calback function and assign to sliders. 
cb = CustomJS(args=args, code=js)
for s in [ref_c_slider, ref_rep_slider, ref_ep_slider, 
          c_slider, rep_slider, ep_slider]:
    s.callback = cb
    s.js_on_change('value', cb)


# Set up things to mimic the slide
header = Div(text="""
<h1 style="font-family:NanumMyeongjo; font-size: 2.5em; 
border-bottom: 1px solid #3c3c3c; width:100%; text-align: left"> Relative changes in parameter values can be mapped to shifts in free energy </h1> 
<br/>
<img src="delF_ref_annotated.png" style="width: 90%; margin-left: 5%; padding-bottom: 10%;">
""")


# Set up the layout
ref_box = widgetbox(ref_c_slider, ref_rep_slider, ref_ep_slider, width=int(plot_width))
box = widgetbox(c_slider, rep_slider, ep_slider, width=int(plot_width))
col = bokeh.layouts.column(ref_box, box)
row1 = bokeh.layouts.row(col, p_fc, p_bohr)
row2 = bokeh.layouts.row(p_dF_c, p_dF_rep, p_dF_ep)
lay = bokeh.layouts.column(header, row1, row2)
bokeh.io.save(lay)



# %%
