### This file contains custom validators used in the project

# Python native imports
import os
import re

# Python 3rd party imports

# Parapy imports
from kbeutils.data import airfoils
from parapy.core.validate import AdaptedValidator

def airfoil_validator(name: str) -> bool:
    naca4_pattern = re.compile(r'^NACA\s?\d{4}$', re.IGNORECASE)
    name_upper = name.upper()

    if naca4_pattern.match(name_upper):
        return True
    else:
        in_input = os.path.isfile(os.path.join(os.getcwd(), "input", "airfoils", name + ".dat"))
        in_kbeutils = os.path.isfile(os.path.join(airfoils.__path__[0], name + ".dat"))
        found = in_input or in_kbeutils
        return found

airfoil_found = AdaptedValidator(airfoil_validator)

def validate_equal_length_lists(object, attribute_names: list[str]):
    list_attributes = [getattr(object, name) for name in attribute_names]
    lengths = [len(lst) for lst in list_attributes]
    lengths_unique = set(lengths)
    if len(lengths_unique) > 1:
        raise ValueError(f"All section-related input lists must have the same length. "
                         f"Got lengths: {lengths} for attributes {attribute_names}")
    return True