### This file contains custom validators used in the project

# Python native imports
import os

# Python 3rd party imports

# Parapy imports
from kbeutils.data import airfoils
from parapy.core.validate import AdaptedValidator

kbeutils_airfoil_found = AdaptedValidator(lambda name: os.path.isfile(os.path.join(airfoils.__path__[0], name + ".dat")))

input_airfoil_found = AdaptedValidator(lambda name: os.path.isfile(os.path.join(os.getcwd(), "input", "airfoils", name + ".dat")))