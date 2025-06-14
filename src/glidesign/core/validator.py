### This file contains custom validators used in the project

# Python native imports
import os
import re

# Python 3rd party imports

# Parapy imports
from kbeutils.data import airfoils
from parapy.core.validate import AdaptedValidator

NACA4_PATTERN = re.compile(r'^NACA\s?(\d{4})$', re.IGNORECASE)
NACA5_PATTERN = re.compile(r'^NACA\s?(\d{5})$', re.IGNORECASE)

def airfoil_validator(name: str) -> bool:
    name_upper = name.upper()

    if NACA4_PATTERN.match(name_upper):
        return True
    elif NACA5_PATTERN.match(name_upper):
        return True
    else:
        in_input = os.path.isfile(os.path.join(os.getcwd(), "input", "airfoils", name + ".dat"))
        in_kbeutils = os.path.isfile(os.path.join(airfoils.__path__[0], name + ".dat"))
        found = in_input or in_kbeutils
        return found

airfoil_found = AdaptedValidator(airfoil_validator)