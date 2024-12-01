#!/usr/bin/env python
import numpy as np


def beta(A: int | float, Z: int | float, Q: float) -> float:
    """Function to calculate deformation parameter `beta` of a nuclei in its ground state.

    Args:
        A (int or float): atomic mass (i.e. total no. nucleons)
        Z (int or float): atomic number (i.e. total number of protons)
        Q (float): quadrupole moment of nucleus in its ground state.  
    """
    R = 1.2 * A**(1/3)
    return np.sqrt(5 * np.pi) / (3 * Z * R**2) * Q

