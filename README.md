# OPFIR modeling

This repository models a THz gas laser for linear molecules such as N2O, HCN, etc.

Structure ParamsLinearMolecule defines the parameters for the problem, with
a constructor function LinearMolecule(). The input for the constructor includes
cavity properties such as the radius, length, temperature, cavity wall material, cavity mode, transmission/reflection coefficients for front and back mirrors. It also includes some physical constants such as Planck constant, speed of light, etc.
Users also need to input gain molecules properties. A simple input of molecular name would choose the correct parameters for the users, but users can also change the parameters if they are not well-known. 
