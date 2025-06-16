import pandas as pd
import os

def load_glider_params(path: str) -> dict:
    """
    Loads glider parameters from an Excel file and stores them in a dictionary.

    Args:
        path (str): The relative path to the Excel file containing the glider parameters.

    Returns:
        dict: A dictionary where the keys are parameter names (as strings) 
              and the values are the corresponding parameter values from the Excel file.
    
    Example:
        params = load_glider_params('glider_params.xlsx')
        print(params)
        # Output: {'span': 10.5, 'chord': 2.0, 'airfoil': 'NACA 0012', ...}
    """
    # Construct the full path to the Excel file using the current working directory
    excel_path = os.path.join(os.getcwd(), path)

    # Read the Excel file into a pandas DataFrame
    df = pd.read_excel(excel_path)

    # Initialize an empty dictionary to store the parameters
    param_dict = {}

    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        key = str(row[0]).strip()  # The first column represents the parameter name
        value = row[1]  # The second column represents the value of the parameter
        param_dict[key] = value  # Add the parameter and its value to the dictionary

    return param_dict