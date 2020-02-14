#%%
import numpy as np
import pandas as pd
import cmdstanpy 

# Load the data
data = pd.read_csv('../../../data/RazoMejia2018_data.csv')
data = data[data['repressors'] > 0]
model = cmdstanpy.CmdStanModel('infer_empirical_F.stan')

#%%

# Define the stan model. 
model_code = """
data { 
    int N; // Number of measurements
    vector fc[N]; // Fold-change values
}

parameters { 
    real<lower=0, upper=1> mu; 
    real<lower=0> sigma;
}

model {
    mu ~ uniform(0, 1);
    sigma ~ normal(0, 0.1);
    fc ~ normal(mu, sigma);
}
"""
model = pystan.StanModel(model_code=model_code)

#%%
# Group by the repressor, operator, and IPTGuM and assign an index


  