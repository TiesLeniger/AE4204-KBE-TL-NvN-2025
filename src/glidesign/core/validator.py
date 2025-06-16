# This file contains custom validators used in the project

# Python native imports
import os
import re

# Python 3rd party imports

# Parapy imports
from kbeutils.data import airfoils
from parapy.core.validate import AdaptedValidator

# Regular expression patterns to match NACA 4-digit and 5-digit airfoil names
NACA4_PATTERN = re.compile(r'^NACA\s?(\d{4})$', re.IGNORECASE)
NACA5_PATTERN = re.compile(r'^NACA\s?(\d{5})$', re.IGNORECASE)

def airfoil_validator(name: str) -> bool:
    """
    Custom validator function to check if a given airfoil name is valid.
    
    The airfoil name is considered valid if it matches the NACA 4-digit or 
    5-digit airfoil naming convention or if a corresponding `.dat` file is 
    found in the 'input/airfoils' directory or within the `kbeutils` library.
    
    Args:
        name (str): The name of the airfoil to be validated.
    
    Returns:
        bool: True if the airfoil name is valid (either a matching NACA code or 
              an existing file), False otherwise.
    
    Example:
        airfoil_validator('NACA 2412')  # True if 'NACA 2412' exists or matches the pattern
    """
    name_upper = name.upper()  # Convert name to uppercase to ensure case-insensitivity

    # Check if the name matches the NACA 4-digit or 5-digit pattern
    if NACA4_PATTERN.match(name_upper):
        return True
    elif NACA5_PATTERN.match(name_upper):
        return True
    else:
        # Check if the airfoil file exists in the input directory or in kbeutils data
        in_input = os.path.isfile(os.path.join(os.getcwd(), "input", "airfoils", name + ".dat"))
        in_kbeutils = os.path.isfile(os.path.join(airfoils.__path__[0], name + ".dat"))
        
        # Return True if the file is found in either location
        found = in_input or in_kbeutils
        return found

# Create a custom validator using AdaptedValidator, which uses the airfoil_validator function
airfoil_found = AdaptedValidator(airfoil_validator)