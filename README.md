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

- install autogluon with pip
````angular2html
pip install autogluon
````
- conda install libs for equibind
````angular2html
conda install pytorch
conda install cudatoolkit=10.2
conda install torchaudio -c pytorch
conda install dgl-cuda10.2 -c dglteam
conda install rdkit -c rdkit
conda install biopandas -c conda-forge
conda install pot -c conda-forge
conda install dgllife -c conda-forge
conda install pyaml -c conda-forge
conda install icecream -c conda-forge
conda install matplotlib -c conda-forge
conda install tensorboard -c conda-forge
````
- install psycopg2
````angular2html
pip install psycopg2
````
or
````angular2html
conda install psycopg2
````
- additional installation may be of help
````angular2html
conda install protobuf=3.20.1
conda install tensorboard=1.15.0
````