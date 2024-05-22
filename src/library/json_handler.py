# imports
    # standard libraries
import json
from typing import Dict
    # user libraries


# functions
def json_writer(file_name: str, dict: Dict[tuple[str, str], float]):
    dict = {str(k): v for k,v in dict.items()}
    with open(file_name, "w") as f:
        json.dump(dict, f, indent=4)

def json_reader(file_name: str):
    with open(file_name, "r") as f:
        return {eval(k): v for k,v in json.load(f).items()}