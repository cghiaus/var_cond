#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 11:11:58 2023

@author: cghiaus

Variable conductivity plane wall
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def fit(θ, λ, deg):
    """
    Parameters
    ----------
    θ : np.array
        Measured temperature, °C.
    λ : np.array
        Measured conductivity, W/(m·K).
    deg : int
        Degree of the fitting polynomial.

    Returns
    -------
    p : np.array
        coefficients of the fitted model λ = a·θ + b or λ = a·θ² + b·θ + c
    R2 : float
        Coefficient of determination: R²
        https://en.m.wikipedia.org/wiki/Coefficient_of_determination
    """
    # Find line of best fit
    p = np.polyfit(θ, λ, deg)

    # Find R² (coefficient of determination)
    λ_predicted = np.polyval(p, θ)
    SStot = ((λ - λ.mean()) ** 2).sum()     # total sum of squares
    SSres = ((λ - λ_predicted) ** 2).sum()  # residuals sum of squares
    R2 = 1 - (SSres / SStot)

    return p, R2


def poly2model(p, θ0, deg):
    """
    Thansforms the parameters of the polynome fitted to data, p,
    λ = a·θ + b or λ = a·θ² + b·θ + c,
    (i.e., λ = p[0]·θ + p[1] or λ = p[0]·θ² + p[1]·θ + p[2]
    into parameters λ0, β, θ0 of the models
    λ = λ0·(1 + β·(θ - θ0)) or λ = λ0·[1 + β·(θ - θ0)²]

    Parameters
    ----------
    p : np.array
        coefficients of the fitted model λ = a·θ + b or λ = a·θ² + b·θ + c
    θ0 : float or NoneType
        Base temperature, °C.
        - `float` if 1st order model.
        - `None` if 2nd order order model
    deg : int
        Degree of the fitted polynomial.

    Returns
    -------
    λ0 : float
        λ for θ = θ₀ in λ = λ₀(1 + β·(θ - θ₀)) or λ = λ₀·[1 + β·(θ - θ₀)²].
    β : float
        Temperature coefficient of thermal conductivity.
    θ0 : float
        Base temperature, °C.

    """
    if deg == 1:
        λ0 = p[1] + p[0] * θ0
        β = p[0] / (p[1] + p[0] * θ0)
    elif deg == 2:
        λ0 = -(p[1]**2 - 4 * p[0] * p[2]) / (4 * p[0])
        θ0 = -p[1] / (2 * p[0])
        β = - 4 * p[0]**2 / (p[1]**2 - 4 * p[0] * p[2])
    else:
        raise ValueError("deg needs to be 1 or 2.")

    return λ0, β, θ0


def fit_plot(θ, λ, p, λ0, β, θ0, deg, material):
    """
    Plot fitted models of 1st or 2nd order to measured data.
    The models are:
        - polynomial λ = a·θ + b or λ = a·θ² + b·θ + c
        - canonical λ₀ = λ₀·(1 + β·(θ - θ₀)) or λ = λ₀·[1 + β·(θ - θ₀)²]

    Parameters
    ----------
    θ : np.array
        Measured temperature, °C.
    λ : np.array
        Measured conductivity, W/(m·K).
    deg : int
        Degree of the fitting polynomial.
    material : str
        Type of material; title of plot.

    Returns
    -------
    None.

    """
    if deg == 1:
        poly = 'λ = a·θ + b'
        model = 'λ = λ₀·(1 + β·(θ - θ₀))'
        λ_model = λ0 * (1 + β * (θ - θ0))
    elif deg == 2:
        poly = 'λ = a·θ² + b·θ + c'
        model = 'λ = λ₀·[1 + β·(θ - θ₀)²]'
        λ_model = λ0 * (1 + β * (θ - θ0)**2)
    else:
        raise ValueError("deg needs to be 1 or 2.")

    λ_poly = np.polyval(p, θ)

    # Plot results
    fig, ax = plt.subplots()
    ax.scatter(θ, λ,
               label='Measured points', marker='o')
    ax.plot(θ, λ_poly,
            label=poly,
            color='blue')
    ax.plot(θ, λ_model,
            label=model,
            color='red')
    ax.set_title(material)
    ax.set_xlabel('Temperature, θ / [°C]')
    ax.set_ylabel('Conductivity, λ / [W/(m·K)]')
    plt.legend()
    plt.show()
    return None


def fit_data2plot(θ0, deg, material):
    """
    Plot data from file and fitted models.

    Parameters
    ----------
    material : TYPE
        DESCRIPTION.
    deg : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    data = '../data/conductivity/'

    df = pd.read_csv(data + material + '.csv')
    θ = df['Temperature / °C']
    λ = df['Conductivity / W/(m·K)']

    p, R2 = fit(θ, λ, deg)
    λ0, β, θ0 = poly2model(p, θ0, deg)
    fit_plot(θ, λ, p, λ0, β, θ0, deg=deg,
             material=material)

    if deg == 1:
        print(f'a = {p[0]:.4f}, b = {p[1]:.4f}')
    elif deg == 2:
        print(f'a = {p[0]:.3e}, b = {p[1]:.3e}, c = {p[2]:.3e}')

    print(f'λ₀ = {λ0:.3f} W/(m·K), β = {β:.5f} W/(m·K²), θ0 = {θ0:.1f} °C')
    print(f'Coefficient of determination: R² = {R2:.4f} ')
    return None


def dirichlet_num(conductivity_model, width, mesh, surf_temp):
    """
    Numerical solution to a plane wall with Dirichlet boundary conditions:
        T0 for x = 0; T1 for x = w
    when conductivity varies with temperature:
        linear λ(T) = λ0·(1 + β·(T - Tb)) if deg == 1
        quadratic λ(T) = λ0·(1 + β·(T - Tb)²) if deg == 2

    Thermal network:
        T0->G[0]->*->G[1]->*...*->G[nb-1]->T1

    Parameters
    ----------
    conductivity_model : list, parameters in λ(T) = λ0·(1 + β·(T - Tb)^deg)
        λ0 : float
            Conductivity for base temperature T = Tb, W/(m·K).
        β : float
            Temperature coefficient in λ = λ(T), W/(m·K²).
        Tb : float
            Base temperature in λ = λ(T), °C.
        deg : int
            Degree of the model: 1 - linear or 2 - quadratic.
        width : float
            Width of the plane wall, m.
        mesh : int
            Number of meshes (i.e., branches) in the numerical model.
    surf_temp : list
        T0 : float
            Temperature on boundary x = 0, °C.
        T1 : float
            Temperature on boundary x = w,, °C.

    Returns
    -------
    θ : np.array
        Temperature at nodes 0, 1, ... , nb-2, °C.
    q : np.array
        Heat flux in branches 0, 1, ... , nb-1, W/m².
    x : np.array
        Locations of temperatures θ on x-axis, m.

    """

    λ0, β, Tb, deg = conductivity_model
    w = width
    nb = mesh
    T0, T1 = surf_temp

    x = np.linspace(0, w, nb)

    A = -np.diff(np.eye(nb))    # incidence matrix
    b = np.zeros(nb)
    b[0], b[-1] = T0, -T1       # temperature source vector
    f = np.zeros(nb - 1)        # flow source vector

    # initial temperature distribution estimated with constant conductivity
    G = λ0 / (w / nb) * np.eye(nb)
    θ0 = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

    ε = 0.001  # °C, acceptable error
    not_convergent = True

    while not_convergent:
        θ_b = np.concatenate(([T0], θ0, [T1]))
        T_0, T_1 = θ_b[:-1], θ_b[1:]    # boundary temp for branches
        if deg == 1:
            θ_mean = (T_0 + T_1) / 2 - Tb
        elif deg == 2:
            θ_mean = (T_0**2 + T_0 * T_1 + T_1**2) / 3 - Tb * (T_0 + T_1 - Tb)
        else:
            raise ValueError("deg needs to be 1 or 2.")

        λ_mean = λ0 * (np.eye(nb) + β * np.diag(θ_mean))
        G = λ_mean / (w / nb)

        θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

        not_convergent = max(abs(θ0 - θ)) > ε
        θ0 = θ

    q = G @ (-A @ θ + b)

    x = np.linspace(0, w, nb + 1)[1: -1]
    return θ, q, x


def dirichlet_anal(cond_model, width, surf_temp, num=50):
    """
    Analytical model

    Parameters
    ----------
    cond_model : list of parameters in λ(T) = λ0·(1 + β·(T - Tb)^deg)
        λ0 : float
            Conductivity for base temperature T = Tb, W/(m·K).
        β : float
            Temperature coefficient in λ = λ(T), W/(m·K²).
        Tb : float
            Base temperature in λ = λ(T), °C.
        deg : int
            Degree of the model: 1 - linear or 2 - quadratic.
    width : float
        Width of the plane wall, m.
    surf_temp : list
        T0 : float
            Temperature on boundary x = 0, °C.
        T1 : float
            Temperature on boundary x = w,, °C.
    num : int
        Number of points in which the model is evaluated.

    Returns
    -------
    θ : np.array
        Temperature at nodes 0, 1, ... , nb-2, °C.
    q : np.array
        Heat flux in branches 0, 1, ... , nb-1, W/m².
    x : np.array
        Locations of temperatures θ on x-axis, m.

    """

    λ0, β, Tb, deg = cond_model
    w = width
    T0, T1 = surf_temp

    x = np.linspace(0, w, num)
    θ0, θ1 = T0 - Tb, T1 - Tb

    if deg == 1:
        C = (θ1 - θ0) / w * (1 + β * ((θ1 + θ0) / 2))
        D = θ0 + β * θ0**2 / 2

        b, c = 1 / β, -2 * (C * x + D) / β
        θ = -b + np.sqrt(b**2 - c)
    elif deg == 2:
        C = (θ1 - θ0) / w * (1 + β * (θ0**2 + θ0 * θ1 + θ1**2) / 3)
        D = θ0 + β * θ0**3 / 3

        b, c = 3 / β, -3 * (C * x + D) / β
        bc = (27 * c / 2 + np.sqrt(108 * b**3 + 729 * c**2) / 2)**(1 / 3)
        θ = b / bc - bc / 3

    θ += Tb
    q = - λ0 * C
    return θ, q, x


def dirichlet_plot(x_n, θ_n, x_a, θ_a, surf_temp, material):
    """

    Plot numerical and analytical solutions of wall with Dirichlet boundary
    conditions.

    Parameters
    ----------
    x_n : np.array
        Locations of numerical solutions θ on x-axis, m.
    θ_n : np.array
        Temperature at nodes in the wall, °C..
    x_a : np.array
        Locations of analytical solutions θ on x-axis, m.
    θ_a : np.array
        Temperatures obtained analytically in the wall, °C.
    surf_temp : list
        T0 : float
            Temperature on boundary x = 0, °C.
        T1 : float
            Temperature on boundary x = w,, °C.
    material : str
        Type of material; title of plot.

    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axes object containing the plotted data.

    """
    T0, T1 = surf_temp

    x_n = np.concatenate(([0], x_n, [x_n[-1] + (x_n[-1] - x_n[-2])]))
    θ_n = np.concatenate(([T0], θ_n, [T1]))

    fig, ax = plt.subplots()
    ax.scatter(x_n, θ_n,
               color='blue', marker='o')
    ax.plot(x_n, θ_n,
            label='Numerical solution', color='blue')
    ax.plot(x_a, θ_a,
            label='Analytical solution', color='orange')
    ax.plot(x_a[[0, -1]], θ_a[[0, -1]],
            label='Linear solution', color='green')
    ax.set_title(material)
    ax.set_xlabel('Wall width, w / [m]')
    ax.set_ylabel('Temperature, $\\theta$ / [°C]')
    plt.legend()
    # plt.show()
    return ax


def wall_num(conductivity_model, width, mesh, air_temp, flow_surf, conv_coef):
    """
    Numerical solution to a plane wall with temperature & flow boundary cond.:
        - temperature T0 of outdoor air; T1 of indoor air
        - flow Q0 for x = 0; Q1 for x = w
        - convection coef. h0 for outdoor air; h1 for indoor air
    when conductivity varies with temperature:
        - linear λ(T) = λ0·(1 + β·(T - Tb)) if deg == 1
        - quadratic λ(T) = λ0·(1 + β·(T - Tb)²) if deg == 2

    Parameters
    ----------
    conductivity_model : list, parameters in λ(T) = λ0·(1 + β·(T - Tb)^deg)
        λ0 : float
            Conductivity for base temperature T = Tb, W/(m·K).
        β : float
            Temperature coefficient in λ = λ(T), W/(m·K²).
        Tb : float
            Base temperature in λ = λ(T), °C.
        deg : int
            Degree of the model: 1 - linear or 2 - quadratic.
        width : float
            Width of the plane wall, m.
        mesh : int
            Number of meshes (i.e., branches) in the numerical model.
    air_temp : list
        T0 : float
            Temperature of outdoor air, °C.
        T1 : float
            Temperature of indoor air, °C.
    flow_surf : list
        Q0 : float
            Flow on surface at x = 0, W/m².
        Q1 : float
            Flow on surface at x = w, W/m².
    conv_coef : list
        h0 : float
            Convection coefficient on surface at x = 0, W/(m²·K).
        h1 : float
            Convection coefficient on surface at x = w, W/(m²·K).

    Returns
    -------
    θ : np.array
        Temperature at nodes 0, 1, ... , nb-2, °C.
    q : np.array
        Heat flux in branches 0, 1, ... , nb-1, W/m².
    x : np.array
        Locations of temperatures θ on x-axis, m.

    """

    λ0, β, Tb, deg = conductivity_model
    w = width
    nb = mesh
    T0, T1 = air_temp
    Q0, Q1 = flow_surf
    h0, h1 = conv_coef

    x = np.linspace(0, w, nb)

    A = -np.diff(np.eye(nb + 2))    # incidence matrix
    b = np.zeros(nb + 2)
    b[0], b[-1] = T0, -T1           # temperature source vector
    f = np.zeros(nb + 1)
    f[0], f[-1] = Q0, Q1            # flow source vector

    # initial temperature distribution estimated with constant conductivity
    Gwall = λ0 / (w / nb) * np.ones(nb)
    G = np.diag(np.concatenate(([h0], Gwall, [h1])))
    θ0 = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

    ε = 0.001  # °C, acceptable error
    not_convergent = True

    while not_convergent:
        θ_b = np.concatenate(([T0], θ0, [T1]))
        T_0, T_1 = θ_b[1:-2], θ_b[2:-1]
        if deg == 1:
            θ_mean = (T_0 + T_1) / 2 - Tb
        elif deg == 2:
            θ_mean = (T_0**2 + T_0 * T_1 + T_1**2) / 3 - Tb * (T_0 + T_1 - Tb)
        else:
            raise ValueError("deg needs to be 1 or 2.")

        λ_mean = λ0 * (1 + β * θ_mean)
        Gwall = λ_mean / (w / nb)
        G = np.diag(np.concatenate(([h0], Gwall, [h1])))
        θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

        not_convergent = max(abs(θ0 - θ)) > ε
        θ0 = θ

    q = G @ (-A @ θ + b)

    x = np.linspace(0, w, nb + 1)
    return θ, q, x


def wall_plot(x, θ, x_a, θ_a, air_temp, material):
    """

    Plot numerical and analytical solutions of plane wall with convection and
    flow sources on the surfaces.

    Parameters
    ----------
    x : np.array
        Locations of numerical solutions θ on x-axis, m.
    θ : np.array
        Temperature at nodes (surfaces and in the meshes of the wall), °C.
    x_a : np.array
        Locations of analytical solutions θ on x-axis, m.
    θ_a : np.array
        Temperatures obtained analitically in the meshes of the wall.
    air_temp : list
        T0 : float
            Temperature outdoor air, °C.
        T1 : float
            Temperature indoor air, °C.
    material : str
        Type of material; title of plot.

    Returns
    -------
    ax : matplotlib.axes._subplots.AxesSubplot
        The axes object containing the plotted data.

    """

    # Plot temperaure in the wall
    x_w, θ_w = x[1:-1], θ[1:-1]
    surf_temp = θ[0], θ[-1]
    ax = dirichlet_plot(x_w, θ_w, x_a, θ_a, surf_temp, material)

    # Plot air temperature
    w_bl = 0.020                            # width of boundary layer

    # location of air temperature
    x_air = np.array(
        [x[0] - 4 * w_bl, x[0] - w_bl,      # outdoor
         x[-1] + w_bl, x[-1] + 4 * w_bl])   # indoor

    T0, T1 = air_temp

    ax.plot(x_air[0:2], [T0, T0],
            label='Air', color='red')
    ax.plot(x_air[2:4], [T1, T1],
            color='red')

    # temperature in boundary layer
    ax.plot([x_air[1], x[0]], [T0, θ[0]],
            color='red', linestyle='dotted')
    ax.plot([x_air[2], x[-1]], [T1, θ[-1]],
            color='red', linestyle='dotted')
    ax.legend()

    ax.axvline(x=x[0], color='gray')
    ax.axvline(x=x[-1], color='gray')
    ax.grid(True)

    return ax


def wall_sim(width, mesh,
             air_temp, flow_surf, conv_coef,
             Tb, conductivity_poly, deg,
             material):
    """
    Simulation for a wall with covection and flow rate on the surfaces.

    Parameters
    ----------
    width : float
        Width of the plane wall, m.
    mesh : int
        Number of meshes (i.e., branches) in the numerical model.
    air_temp : list
        T0 : float
            Temperature of outdoor air, °C.
        T1 : float
            Temperature of indoor air, °C.
    flow_surf : list
        Q0 : float
            Flow on surface at x = 0, W/m².
        Q1 : float
            Flow on surface at x = w, W/m².
    conv_coef : list
        h0 : float
            Convection coefficient on surface at x = 0, W/(m²·K).
        h1 : float
            Convection coefficient on surface at x = w, W/(m²·K).
    Tb : float
        Base temperature in λ = λ(T), °C.
    conductivity_poly : list of float
        - [a, b] if fitted polynom is linear, λ(T) = a·T + b.
        - [a, b, c] if titted polynom is quadratic, λ(T) = a·T² + b·T + c.
    deg : int
        - 1 for linear model, λ(T) = a·T + b.
        - 2 for quadratic model, λ(T) = a·T² + b·T + c.
    material : string
        Name of the material.

    Returns
    -------
    q : np.array
        Heat flux in branches 0, 1, ... , nb-1, W/m².
    conductivity_model : list
        - λ0, β, Tb: λ(T) = λ0·(1 + β·(T - Tb)) or λ(T) = λ0·(1 + β·(T - Tb)²)
        - deg is 1 (linear) or 2 (quadratic)
    surf_temp : list of float
        - θ[0] temperature at wall surface x = 0, °C.
        - θ[-1] temperature at wall surface x = w, °C.

    """

    # λ(T) = a·T + b to λ(T) = λ0·(1 + β·(T - Tb))
    # or
    # λ(T) = a·T² + b·T + c to λ(T) = λ0·(1 + β·(T - Tb)²)
    λ0, β, Tb = poly2model(conductivity_poly, Tb, deg)
    conductivity_model = [λ0, β, Tb, deg]

    θ, q, x = wall_num(conductivity_model, width, mesh,
                       air_temp, flow_surf, conv_coef)
    # x_w, θ_w = x[1:-1], θ[1:-1]     # wall

    surf_temp = θ[0], θ[-1]

    # Analytical solution
    θ_a, q_a, x_a = dirichlet_anal(conductivity_model, width,
                                   surf_temp, num=5)

    wall_plot(x, θ, x_a, θ_a, air_temp, material)

    return q, conductivity_model, surf_temp
