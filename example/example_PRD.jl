"""
An example script to plot a payload-range diagram with a 
fleet of missions
"""

# 1. Import modules
using TASOPT
using Plots
# import indices for calling parameters

# Load default model
ac = load_default_model() #Use default model for payload-range diagram
size_aircraft!(ac)

#Create Payload-Range diagram
TASOPT.PayloadRange(ac) #This function is in src/outputs.jl