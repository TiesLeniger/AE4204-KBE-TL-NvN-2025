import numpy as np

def convert_matlab_dict(matlab_dict: dict) -> dict:
    """
    Converts a MATLAB-style dictionary (with numpy arrays or lists as values) 
    to a Python dictionary where each value is flattened into a list.
    
    Args:
        matlab_dict (dict): A dictionary where each key corresponds to a MATLAB 
                             variable, and each value is typically a numpy array 
                             or a list that may have multiple dimensions.
    
    Returns:
        dict: A new dictionary with the same keys, but with each value converted 
              to a flattened list.
    
    Example:
        matlab_dict = {'a': np.array([[1, 2], [3, 4]]), 'b': [5, 6, 7]}
        result = convert_matlab_dict(matlab_dict)
        print(result)
        # Output: {'a': [1, 2, 3, 4], 'b': [5, 6, 7]}
    """
    # Create a new dictionary by iterating over the input dictionary
    return {
        key: np.array(value).flatten().tolist()  # Flatten the numpy array and convert to list
        for key, value in matlab_dict.items()  # Iterate through each key-value pair in the input dictionary
    }