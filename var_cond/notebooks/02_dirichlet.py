#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 18:03:14 2023

@author: cghiaus
"""

import sys
sys.path.append('../src')  # Add 'scr' directory to the system path
import var_cond

# Data
w = 0.20                    # m, width of the plane wall
nb = 3                      # number of branches, nb > 1
T0, T1 = -20, 50            # °C, boundary temperatures out/in surface

