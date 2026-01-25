import os
import re
import json
from itertools import product

import os
import argparse

def makedir_if_not_exists(directory):
    os.makedirs(directory, exist_ok=True)

def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        print(f"{dirname} is not a valid directory.")
        exit
    else:
        return dirname

def is_file(filename):
    """Checks if the provided file exists and is accessible."""
    if not os.path.isfile(filename):  # Check if it's a file
        print(f"{filename} is not a valid file.")
        exit
    return filename  # Return the valid filename as the value of the argument