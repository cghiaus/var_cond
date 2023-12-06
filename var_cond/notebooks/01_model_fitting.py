#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 20:13:08 2023

@author: cghiaus

Folder structure:
var_cond
    |--data
        |--conductivity
            |--Fiberglass.csv
            |--Polyisocyanurate.csv
    |--notebooks
        |--01.ipynb
    |--src
        |--__iniy__.py
        |--var_cond.py

"""

import pandas as pd
import sys
sys.path.append('../src')  # Adding the 'scr' directory to the system path
import var_cond

data = '../data/conductivity/'

material = 'Fiberglass'
df = pd.read_csv(data + material + '.csv')

θ = df['Temperature / °C']
λ = df['Conductivity / W/(m·K)']
θ0 = 10     # °C, base temperature
deg = 1     # model degree

[a, b], R2 = var_cond.fit(θ, λ, deg)
λ0, β, θ0 = var_cond.poly2model([a, b], θ0, deg)
var_cond.fit_plot(θ, λ, [a, b], λ0, β, θ0, deg=deg,
                  material=material)

material = 'Polyisocyanurate'
df = pd.read_csv(data + material + '.csv')

θ = df['Temperature / °C']
λ = df['Conductivity / W/(m·K)']
deg = 2     # model degree

[a, b, c], R2 = var_cond.fit(θ, λ, deg)
λ0, β, θ0 = var_cond.poly2model([a, b, c], θ0, deg)
var_cond.fit_plot(θ, λ, [a, b, c], λ0, β, θ0, deg=deg,
                  material=material)
