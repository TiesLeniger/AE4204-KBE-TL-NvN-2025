import numpy as np

def convert_matlab_dict(matlab_dict: dict) -> dict:
    return {
        key: np.array(value).flatten().tolist()
        for key, value in matlab_dict.items()
    }