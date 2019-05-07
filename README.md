# SemiLagPy - Semi Lagrangian scheme in python

This repository shows how to implement a Semi-Lagrangian scheme in python.
During a course in Numerical Weather Prediction (NWP), I found it quiet hard to find some straightforward examples on how to implement such a scheme. The scheme is widely used in NWP and can be illustrated using a simple 1D advection model. The major advantage of this scheme free choice in the time step. This step size can be bigger than in normal finite difference schemes (e.g. Leap-Frog).

The [src/SemiLagrangian.py](https://github.com/christian512/SemiLagPy/blob/master/src/SemiLagrangian.py) already includes a prototype class for integrating periodic, one dimensional problems using a linear or cubic interpolation.
An example of using this class is given in [advection1D.ipynb](https://github.com/christian512/SemiLagPy/blob/master/advection1D.ipynb)
