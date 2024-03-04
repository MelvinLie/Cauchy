# Cauchy
A lightweight boundary element method (Cauchy method) in python and C.

This code was developed for the standard analysis of single stretched wire measurements.

See the publication 'Boundary-Element Methods for Field Reconstruction in Accelerator Magnets' by
M. Liebsch for details.

[Link to the publication](https://ieeexplore.ieee.org/document/8954937)

This code is using the poetry packaging manager. All You need to do is to run:
1. poetry build
2. poetry install

Then You can launch the examples by

poertry run examples/plot_example_1.py

I advice You to generate the autodocs. To do so, navigate to the docs directory and launch

poetry run make html

Autor: Melvin Liebsch

Emai: melvin.liebsch@cern.ch
