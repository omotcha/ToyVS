"""
platform: win
env: any
name: config_win.py
configurations of ToyVS
"""

import os
project_dir = os.getcwd()
data_dir = os.path.join(project_dir, "data")
# in this configuration, 2 csv files with nearly 2k and 8m lines of smiles strings are used as data source
data_2k = os.path.join(data_dir, "smiles_2k.csv")
data_8m = os.path.join(data_dir, "smiles_8m.csv")


if __name__ == '__main__':
    pass
