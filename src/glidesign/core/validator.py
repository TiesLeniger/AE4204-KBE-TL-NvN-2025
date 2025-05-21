### This file contains custom validators used in the project

# Python native imports
import os

# Python 3rd party imports

# Parapy imports
from kbeutils.data import airfoils
from parapy.core.validate import AdaptedValidator

airfoil_found = AdaptedValidator(lambda name: (
        os.path.isfile(os.path.join(os.getcwd(), "input", "airfoils", name + ".dat")) or
        os.path.isfile(os.path.join(airfoils.__path__[0], name + ".dat"))
    ))

def validate_equal_length_lists(object, attribute_names: list[str]):
    list_attributes = [getattr(object, name) for name in attribute_names]
    lengths = [len(lst) for lst in list_attributes]
    lengths_unique = set(lengths)
    if len(lengths_unique) > 1:
        raise ValueError(f"All section-related input lists must have the same length. "
                         f"Got lengths: {lengths} for attributes {attribute_names}")
    return True