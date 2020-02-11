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
bokeh.plotting.output_file("../assets/bohr_explorer.html")

# Define the various sliders. 
rep_slider = Slider(title='repressors per cell', start=0, end=1000, step=10, 
                    value=200, bar_color=colors['light_red'])
ep_slider = Slider(title='DNA binding energy [kT]', start=-25, end=0, 
                    step=0.5, value=-14, bar_color=colors['light_orange'])
c_slider = Slider(title='log\u2081\u2080 inducer concentration', start=-4, end=6, step=0.01, value=0, bar_color=colors['light_blue'])


# Define the range of inducer concentrations
c_range = np.logspace(-4, 6, 200)

# Compute the starting position for point and curve. 
fc = phd.thermo.SimpleRepression(R=rep_slider.value, ep_r=ep_slider.value,
           ka=200, ki=1, ep_ai=5, effector_conc=c_range).fold_change()
fc_point = phd.thermo.SimpleRepression(R=rep_slider.value, ep_r=ep_slider.value,
           ka=200, ki=1, ep_ai=5, effector_conc=10**c_slider.value).fold_change()
bohr_point = phd.thermo.SimpleRepression(R=rep_slider.value, ep_r=ep_slider.value,
           ka=200, ki=1, ep_ai=5, effector_conc=10**c_slider.value).bohr_parameter()

# Set up the data source. 
source = ColumnDataSource(pd.DataFrame({'fc':fc, 'c':c_range, 'c_ki':np.log10(c_range)}))
point_source = ColumnDataSource({'c':[10**c_slider.value], 'c_ki':[c_slider.value],
                                 'fc':[fc_point], 'bohr': [bohr_point]})


# Define the callback to update the sources. 
js = """
var rep = rep_slider.value;
var c = Math.pow(10, c_slider.value);
var ep_r = ep_slider.value;
var curve = fc_source.data;
var point = point_source.data;
var c_range = curve['c'];
var ka = ka;
var ki = ki;
var ep_ai = ep_ai;
var n = 2;


// Define a function to simplify caluclation of things. 
function computePact(c, ka, ki, ep_ai, n) {
        var numer = Math.pow(1 + c / ka, n);
        var denom = numer + Math.exp(-ep_ai) * Math.pow(1 + c/ ki, n);
        return numer / denom
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


// Given the slider data, update the curve and point
var updated_point_fc = computeFoldChange(rep, ep_r, c, ka, ki, ep_ai, n);
var updated_point_bohr = computeBohrParameter(rep, ep_r, c, ka, ki, ep_ai, n)
point_source.data['fc'] = [updated_point_fc];
point_source.data['bohr'] = [updated_point_bohr];
point_source.data['c'] = [c];
point_source.data['c_ki'] = [Math.log10(c / ki)];

// Loop through the c range and update the curve. 
for (var i = 0; i < fc_source.data['c'].length; i++) { 
fc_source.data['fc'][i] = computeFoldChange(rep, ep_r, fc_source.data['c'][i], ka, ki, ep_ai, n);
}
// Emit the changes. 
fc_source.change.emit();
point_source.change.emit();
"""

# Assemble the arguments 
args = {'rep_slider':rep_slider, 'c_slider':c_slider, 'ep_slider':ep_slider, 
        'fc_source':source, 'point_source':point_source, 'ka':200, 'ki':1,
        'ep_ai':5}

# Define the callback and assign to sliders. 
cb = CustomJS(args=args, code=js)

for s in [rep_slider, c_slider, ep_slider]:
    s.callback = cb
    s.js_on_change('value', cb)


# Define the master curve. 
bohr_range = np.linspace(-10, 10, 200)
master_curve = (1 + np.exp(-bohr_range))**-1


# Set up the figure canvases. 
p_fc = bokeh.plotting.figure(width=300, height=300,
                            x_axis_label='log\u2081\u2080(c/ Ki)',
                            y_axis_label='fold-change',
                            y_range=[-0.05, 1.1],
                            title='phenotypic output')
p_bohr = bokeh.plotting.figure(width=300,
                            height=300, x_axis_label='free energy [kT]',
                            y_axis_label='fold-change',
                            title='energetic representation')

# Set the correct font sizes. 
for p in [p_fc, p_bohr]:
    p.axis.axis_label_text_font_size = "1.5em"
    p.title.text_font_size = "1.5em"
    p.axis.major_label_text_font_size = "1em"

p_fc.line(x='c_ki', y='fc', source=source, line_width = 2, color=colors['light_purple'])
p_fc.circle(x='c_ki', y='fc', source=point_source, line_width=2, color=colors['light_purple'], size=8)

p_bohr.line(x=bohr_range, y=master_curve, color=colors['black'], line_width=2)
p_bohr.circle(x='bohr', y='fc', source=point_source, line_width=2, color=colors['purple'], size=8)

# Set up the widget box and layout. 
box = widgetbox(rep_slider, ep_slider, c_slider)
plot_row= bokeh.layouts.row(p_fc, p_bohr, sizing_mode='scale_width')
slider_row = bokeh.layouts.row(c_slider, rep_slider, ep_slider, sizing_mode='scale_width')
lay = bokeh.layouts.column(slider_row, plot_row, sizing_mode='scale_width')
bokeh.io.save(lay)
# %%
