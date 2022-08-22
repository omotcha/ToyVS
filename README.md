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

- install libs for equibind
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

- install postgresql and psycopg2
````angular2html
conda install rdkit-postgresql -c rdkit
pip install psycopg2
````
- install jupyter and nglview for visualisation
````angular2html
conda install jupyter -c conda-forge
conda install nglview -c conda-forge
````

### Requirements for ToyVS + Autogluon(For Test)
In my test, Autogluon generates an LightGBXMT model that predicts faster and better than GBT.
Be aware that Autogluon runs on cuda11.6. 
Installing autogluon directly on previous conda environment may cause conflicts
- create conda env: 
````angular2html
conda create -n toyvs_ag python=3.7
conda activate toyvs_ag
````
- install libs for autogluon
````angular2html
conda install pytorch torchvision torchaudio cudatoolkit=11.6 -c pytorch -c conda-forge
pip install autogluon
````
- install libs for equibind
````angular2html
conda install dgl-cuda11.6 -c dglteam
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
- install postgresql and psycopg2
````angular2html
conda install rdkit-postgresql -c rdkit
pip install psycopg2
````
- install jupyter and nglview for visualisation
````angular2html
conda install jupyter -c conda-forge
conda install nglview -c conda-forge
````

### Database Preparation
- you might want to initialize database
````angular2html
initdb -D [database dir]
pg_ctl -D [database dir] -l logfile start
createdb toydb
psql toydb
````
example may be helpful to you: db create user and grant privileges
````angular2html
CREATE USER omotcha WITH PASSWORD '123456'  VALID UNTIL '2024-01-01';
GRANT CONNECT ON DATABASE toydb TO omotcha;
GRANT USAGE ON SCHEMA public TO omotcha;
GRANT SELECT,UPDATE,INSERT,DELETE ON ALL TABLES IN SCHEMA public TO omotcha;
ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO omotcha;
````

- you might want to load data file like smiles.csv to database
- you can check private member functions implemented in dbutil/dbutil.py and run them in dbtest, like following:
````angular2html
DBUtil._create_table_smiles2k()
DBUtil._create_table_results2k()
DBUtil._fill_table_with_local_file(data2k, 'smiles2k')
DBUtil._create_distinct_table('smiles2k')
````
