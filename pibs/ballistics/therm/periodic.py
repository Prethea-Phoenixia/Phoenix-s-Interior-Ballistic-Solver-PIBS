import json
import os

with open(os.path.join(os.path.dirname(__file__), "periodic.json"), "r") as f:
    molarMasses = json.load(f)
