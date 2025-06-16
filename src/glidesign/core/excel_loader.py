import pandas as pd
from pathlib import Path

def load_glider_params(filename: str) -> dict:
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
    if not Path(filename).suffix == ".xlsx":
        filename += ".xlsx"
    # Construct the full path to the Excel file using the current working directory
    excel_path = Path.cwd() / "input" / filename

    # Read the Excel file into a pandas DataFrame
    df = pd.read_excel(excel_path, sheet_name="Input")

    parameter_list = df["Parameter"].to_list()
    value_list = df["Value"].to_list()

    param_dict = {key: value for key, value in zip(parameter_list, value_list)}

    return param_dict