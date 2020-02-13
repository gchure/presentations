# -*- coding: utf-8 -*-
#%%
import numpy as np
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
bokeh.plotting.output_file("../figs/fixed_wt_model_explorer.html", mode='inline')


# Define parameter ranges
c_range = np.logspace(-3, 6, 300)
subsamp = list(np.arange(0, len(c_range), 20))
bohr_range = np.linspace(-20, 20, 500)

ref_R = 300 
ref_epRA = -14
ref_Ka = 200
ref_Ki = 1
ref_epAI = 5

# Set the reference induction profile
ref_arch = phd.thermo.SimpleRepression(R=ref_R, ep_r=ref_epRA, effector_conc=c_range,
                                    ka=ref_Ka, ki=ref_Ki, ep_ai=ref_epAI, 
                                    n_sites=2)
ref_fc = ref_arch.fold_change()
ref_bohr = ref_arch.bohr_parameter()
ref_delta_bohr = ref_bohr - ref_bohr

# Define the source
source = ColumnDataSource(data=dict(c=c_range, c_Ka=np.log(c_range/ref_Ka), ref_fc=ref_fc, 
                                        mut_fc=ref_fc, mut_bohr=ref_bohr,
                                        mut_delta_bohr=ref_delta_bohr, ref_bohr=ref_bohr))
view = CDSView(source=source, filters=[IndexFilter(subsamp)])

# Instantiate the figure canvas
p_fc = bokeh.plotting.figure(width=500, height=300,
                            x_axis_label='log [c / Ka (reference)]', y_axis_label='fold-change',
                            title='induction profile', x_range=[-6, 5],
                            y_range=[-0.1, 1.1])

p_bohr = bokeh.plotting.figure(width=500, height=300,
                                x_axis_label='free energy [kT]', 
                                y_axis_label='fold-change',
                                x_range=[-15, 15], y_range=[-0.1, 1.1],
                                title='energetic representation')
p_delBohr = bokeh.plotting.figure(width=500, height=300,
                                  x_axis_label='log [c / Ka (reference)]', y_axis_label='∆F [kT]',
                                  y_range=[-15, 15],
                                  title='shift in free energy')

# Plot fold-change values
p_fc.line(x='c_Ka', y='ref_fc', source=source, color=colors['black'], line_width=2)
p_fc.line(x='c_Ka', y='mut_fc', source=source, color=colors['green'], line_width=2)
p_fc.circle(x='c_Ka', y='mut_fc', source=source, view=view,
            color='white', line_width=2, size=10, 
            line_color=colors['green'])

# Plot the data collapse
p_bohr.line(bohr_range, (1 + np.exp(-bohr_range))**-1, color='black', 
            line_width=2)
p_bohr.circle(x='mut_bohr', y='mut_fc', source=source, view=view,
            color='white', line_color=colors['green'], size=10, 
            line_width=2)



# Pot the delta bohr
dbohr_ref = bokeh.models.Span(location=0, dimension='width', line_color='black',
                               line_width=1)
p_delBohr.line(x='c_Ka', y='mut_delta_bohr', source=source, color=colors['green'], 
            line_width=2)
p_delBohr.circle(x='c_Ka' , y='mut_delta_bohr', view=view, source=source,
                color='white',  size=10, 
                line_color=colors['green'], line_width=2)
dbohr_ref = bokeh.models.Span(location=0, dimension='width', line_color='black',
                               line_width=1)
p_delBohr.renderers.extend([dbohr_ref])


# Position legends
# p_bohr.legend.location = 'top_left'
# p_delBohr.legend.location = 'top_left'
# p_fc.legend.location = 'top_left'


# #######################
# CONTROLS
# #######################
# Reference control

# Mutant controls
mut_epRA_slider = Slider(start=-10, end=10, value=0, step=0.1,
title='difference in DNA binding energy [kT]', 
         bar_color=colors['light_orange'])
mut_R_slider = Slider(start=-299, end=1500, value=260, step=1,
 title='difference in repressor expression [rep. per cell]', 
        bar_color=colors['light_red'])
mut_ka_slider = Slider(start=np.log(0.1), end=4, value=0, step=0.01,
    title='log relative change in Ka', bar_color=colors['light_purple'])

mut_ki_slider = Slider(start=-3, end=3, value=0, step=0.01,
    title='log relative change in Ki', bar_color=colors['light_purple'])
            
mut_epAI_slider = Slider(start=-5, end=5, value=0, step=0.1,
    title='difference in allosteric state energy [kT]', 
    bar_color=colors['light_purple'])

    
callback_args = {'source':source,
                 'refepRA':ref_epRA,
                 'refR': ref_R,
                 'refKa': ref_Ka,
                 'refKi': ref_Ki,
                 'refepAI':ref_epAI,
                 'mutepRA':mut_epRA_slider,
                 'mutR': mut_R_slider,
                 'mutKa': mut_ka_slider,
                 'mutKi': mut_ki_slider,
                 'mutepAI':mut_epAI_slider}

reset_mut = CustomJS(args=callback_args, code="""
        mutepRA.value = 0;
        mutR.value = 0;
        mutKa.value = 0;
        mutKi.value = 0;
        mutepAI.value = 0;
        """)

callback  = CustomJS(args=callback_args, code="""
                // Define constants
                var data = source.data;
                var c = data['c'];

                // Reference parameters
                var ref_ep_r = refepRA;
                var ref_R = refR;
                var ref_ka = refKa;
                var ref_ki = refKi;
                var ref_ep_ai= refepAI;
                var ref_n = 2;

                // Mktant parameters
                var mut_fc = data['mut_fc'];
                var mut_bohr = data['mut_bohr'];
                var mut_delta_bohr = data['mut_delta_bohr'];
                var mut_ep_r = mutepRA.value + refepRA;
                var mut_R = mutR.value + refR;
                var mut_ka = Math.exp(mutKa.value) * refKa;
                var mut_ki = Math.exp(mutKi.value) * refKi;
                var mut_ep_ai = mutepAI.value + refepAI;
                var mut_n = 2;

                // Define functions for calculating various quantities
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

                // Evaluate the fold-change of the reference and perturbed state
                for (var i = 0; i < c.length; i++) {
                    mut_fc[i] = computeFoldChange(mut_R, mut_ep_r, c[i], mut_ka, 
                                                  mut_ki, mut_ep_ai, mut_n);
                    mut_bohr[i] =  computeBohrParameter(mut_R, mut_ep_r, c[i],
                                                        mut_ka, mut_ki, mut_ep_ai,
                                                        mut_n);
                    var ref_bohr = computeBohrParameter(ref_R, ref_ep_r, c[i],
                                                        ref_ka, ref_ki, ref_ep_ai,
                                                        ref_n);
                    mut_delta_bohr[i] = mut_bohr[i] - ref_bohr;
                }
                source.change.emit();
                """)

# Define the buttons
mut_reset = Button(label='double-click to reset set mutant to wild type', callback=reset_mut)
mut_reset.js_on_click(callback)

# Assemble controls
mut_controls = [mut_reset, mut_R_slider,  mut_epRA_slider, mut_epAI_slider, mut_ka_slider, mut_ki_slider,
                mut_epAI_slider]

for  mc in mut_controls[1:]:
    mc.callback = callback


header = Div(text="""
<h1 style="font-family:NanumMyeongjo; width: 1000px; font-size:2.5em;border-bottom: 1px solid #3c3c3c; text-align: left;">
The profile of the free energy difference is parameter dependent</h1>
<br/>
<center>
<img src="delF_annotated.png" style="margin-left: 0%; width: 600px; padding-bottom: 15%;"
</center>

"""
)

box = widgetbox(mut_reset, mut_ka_slider, mut_ki_slider, mut_epAI_slider, mut_R_slider, mut_epRA_slider,
width=500, height=300)
row1 = bokeh.layouts.row(box, p_fc)
row2 = bokeh.layouts.row(p_bohr, p_delBohr)
layout = bokeh.layouts.column(header, row1, row2)

                               
bokeh.io.save(layout)



# %%
