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
createdb toydb
psql toydb
````
db create user and grant privileges
````angular2html
CREATE USER omotcha WITH PASSWORD '123456'  VALID UNTIL '2024-01-01';
GRANT CONNECT ON DATABASE toydb TO omotcha;
GRANT USAGE ON SCHEMA public TO omotcha;
GRANT SELECT,UPDATE,INSERT,DELETE ON ALL TABLES IN SCHEMA public TO omotcha;
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO omotcha;
````

- you might want to load data file like smiles.csv to database
- you can check private member functions implemented in dbutil/dbutil.py and run them in dbtest
````angular2html
DBUtil._create_table_smiles2k()
DBUtil._create_table_results2k()
DBUtil._fill_table_with_local_file(data2k, 'smiles2k')
DBUtil._create_distinct_table('smiles2k')
````
