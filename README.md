# ToyVS
A demo VS framework that incorporates EquiBind and ECIF. Slowly in progress...





### References
- EquiBind
https://github.com/HannesStark/EquiBind


- ECIF
https://github.com/DIFACQUIM/ECIF

### Requirements
- create conda env: 
````angular2html
conda create -n toyvs python=3.7
conda activate toyvs
````

- conda install libs for equibind
````angular2html
conda install pytorch
conda install cudatoolkit=10.2
conda install torchaudio -c pytorch
conda install dgl-cuda10.2 -c dglteam
conda install rdkit -c conda-forge
conda install openbabel -c conda-forge
conda install biopython -c conda-forge
conda install biopandas -c conda-forge
conda install pot -c conda-forge
conda install dgllife -c conda-forge
conda install pyaml -c conda-forge
conda install icecream -c conda-forge
conda install matplotlib -c conda-forge
conda install tensorboard -c conda-forge
````

- install postgresql and pscopg2
````angular2html
conda install rdkit-postgresql -c rdkit
pip install psycopg2
````

- you might want to initialize database
````angular2html
initdb -D [database dir]
pg_ctl -D [database dir] -l logfile start
````

- you might want to load data file like smiles.csv to database
- you can check private member functions implemented in dbutil/dbutil.py
````angular2html
DBUtil._create_table_smiles2k()
DBUtil._create_table_smiles8m()
DBUtil._fill_table_with_local_file()
````
