import pandas as pd
import os

def load_glider_params(path: str) -> dict:
    excel_path = os.path.join(os.getcwd(), path)
    df = pd.read_excel(excel_path)

    # Assumes two columns: 'Parameter' and 'Value'
    param_dict = {}
    for _, row in df.iterrows():
        key = str(row[0]).strip()
        value = row[1]
        param_dict[key] = value
    return param_dict