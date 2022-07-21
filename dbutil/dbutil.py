"""
platform: win
env: any but with psycopg2 lib
name: dbutil.py
database: you have to create it first locally
enable connection to db and some basic operations
"""
import os.path
from configs.config_win import *
import pickle
import psycopg2
from rdkit import Chem


class DBUtil:

    def _get_all_tables(self):
        """
        get all the tables in db
        :return:
        """
        query = '''
        SELECT tablename
        FROM pg_tables
        WHERE tablename NOT LIKE 'pg%'
                        AND tablename NOT LIKE 'sql_%'
        ORDER BY tablename; 
        '''
        self._cursor.execute(query)
        result = self._cursor.fetchall()
        return result

    def _connectdb(self):
        """
        connect to local database
        :return:
        """
        self._connection = psycopg2.connect(database='toydb',
                                            user='omotcha',
                                            password='123456',
                                            port='5432',
                                            host='127.0.0.1')
        self._cursor = self._connection.cursor()
        self._tables = [i[0] for i in self._get_all_tables()]

    def __init__(self):
        self._connectdb()

    def __del__(self):
        self._connection.close()

    # private members
    _connection = None
    _cursor = None
    _tables = None

    # private util functions
    def _create_table_smiles2k(self):
        """
        create a table containing 2k smiles strings
        for test
        private cuz only done once
        :return:
        """
        query = '''
        CREATE TABLE smiles2k(
        id serial primary key,
        smiles text,
        molwt float); 
        '''
        self._cursor.execute(query)
        self._connection.commit()
        self._tables = [i[0] for i in self._get_all_tables()]

    def _create_table_results2k(self):
        """
        create a table containing 2k results
        for test
        private cuz only done once
        :return:
        """
        query = '''
        CREATE TABLE results2k(
        id serial primary key,
        equimol BYTEA,
        prediction float); 
        '''
        self._cursor.execute(query)
        self._connection.commit()
        self._tables = [i[0] for i in self._get_all_tables()]

    def _create_table_results8m(self):
        """
        create a table containing 8m results
        for test
        private cuz only done once
        :return:
        """
        query = '''
        CREATE TABLE results8m(
        id serial primary key,
        equimol BYTEA,
        prediction float); 
        '''
        self._cursor.execute(query)
        self._connection.commit()
        self._tables = [i[0] for i in self._get_all_tables()]

    def _create_table_smiles8m(self):
        """
        create a table containing 8m smiles strings
        for test
        private cuz only done once
        :return:
        """
        query = '''
        CREATE TABLE smiles8m(
        id serial primary key,
        smiles text,
        molwt float); 
        '''
        self._cursor.execute(query)
        self._connection.commit()
        self._tables = [i[0] for i in self._get_all_tables()]

    def _drop_table(self, table_name):
        """
        drop table by name
        for test
        :param table_name: the name of table to be dropped
        :return:
        """
        if table_name not in self._tables:
            print('ERROR: table name not found in database')
            return

        query = '''
        DROP TABLE {};
        '''.format(table_name)
        self._cursor.execute(query)
        self._connection.commit()
        self._tables = [i[0] for i in self._get_all_tables()]

    def _fill_table_with_local_file(self, f=data_2k, table_name='smiles2k'):
        """
        fill the table by name with local file by dir
        for test
        private cuz only done once for each local file
        :param f: file dir add file with .csv extension
        :param table_name: the name of table to be dropped
        :return:
        """
        if not os.path.isfile(f):
            print('ERROR: source file not found in local')
            return

        if table_name not in self._tables:
            print('ERROR: table name not found in database')
            return

        with open(f, 'r') as ff:
            next(ff)
            self._cursor.copy_from(ff, table_name, sep=',', )
        self._connection.commit()

    def _create_initial_mol(self, table_name='smiles2k'):
        """
        generate initial molecule table from smiles table
        private cuz only done once for each table
        :param table_name: the name of table to be generated from
        :return:
        """
        if table_name not in self._tables:
            print('ERROR: table name not found in database')
            return

        query = '''
        SELECT * INTO mols_{} FROM (SELECT id, mol_from_smiles(smiles::cstring) m FROM {}) tmp WHERE m IS NOT NULL;
        '''.format(table_name, table_name)
        self._cursor.execute(query)
        self._connection.commit()

    def _test_mol_table(self):
        """
        a test for inserting a
        :return:
        """
        smiles = "n1(c(nnn1)SCC(=O)Nc1sc(nn1)SCC(=O)NNC(=O)CSc1n(nnn1)c1ccccc1)c1ccccc1"
        self._cursor.execute("SELECT equimol FROM results2k WHERE id='{}'".format(1))
        result = pickle.loads(self._cursor.fetchone()[0])
        origin = Chem.MolFromSmiles(smiles)
        block_optimized = Chem.MolToMolBlock(origin)
        with open(os.path.join(output_dir, '1.sdf'), "w") as newfile:
            newfile.write(block_optimized)
        block_optimized = Chem.MolToMolBlock(result)
        with open(os.path.join(output_dir, '2.sdf'), "w") as newfile:
            newfile.write(block_optimized)

    # callables
    def get_cursor(self):
        """
        get the cursor
        :return: cursor
        """
        return self._cursor

    def get_num_rows(self, table_name):
        """
        get the number of rows in table
        :param table_name: the name of table to be queried
        :return:
        """
        if table_name not in self._tables:
            print('ERROR: table name not found in database')
            return

        query = '''
        SELECT count(*) FROM {}
        '''.format(table_name)
        self._cursor.execute(query)
        result = self._cursor.fetchall()
        return result[0][0]

    def fetch_smiles_by_index(self, i, table_name='smiles2k'):
        """
        fetch smiles from smiles table by index
        :param i: index
        :param table_name: smiles table to be queried
        :return: smiles
        """
        if table_name not in self._tables:
            print('ERROR: table name not found in database')
            return

        query = '''
        SELECT smiles FROM {} WHERE id = {}
        '''.format(table_name, i)
        self._cursor.execute(query)
        result = self._cursor.fetchall()[0][0]
        return result

    def dbtest(self):
        """
        for personal tests
        :return:
        """
        # self._drop_table('mols_smiles2k')
        # self._create_initial_mol('smiles2k')
        # self.fetch_mol_by_index(0, 'mols_smiles2k')
        # self._drop_table('test')
        self._drop_table('results2k')
        self._create_table_results2k()

    def insert_prediction(self, i, mol, pred, table_name='results2k'):
        """
        insert one row of result into table
        :param i: index
        :param mol: rdkit molecule object to be pickled
        :param pred: the prediction
        :param table_name: name of table to be inserted
        :return:
        """
        if table_name not in self._tables:
            print('ERROR: table name not found in database')
            return
        sql = "INSERT INTO {}".format(table_name)
        sql = sql + " VALUES(%s, %s, %s)"
        self._cursor.execute(sql, (i, pickle.dumps(mol), pred))
        self._connection.commit()

    def reset_results2k(self):
        """

        :return:
        """
        self._drop_table('results2k')
        self._create_table_results2k()


if __name__ == '__main__':
    dbu = DBUtil()
    dbu.reset_results2k()
