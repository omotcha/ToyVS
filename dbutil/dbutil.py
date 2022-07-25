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

trans_map = {
   '\a': r'\a',
   '\b': r'\b',
   '\f': r'\f',
   '\n': r'\n',
   '\r': r'\r',
   '\t': r'\t',
   '\v': r'\v',
   '\'': r'\'',
   '\"': r'\"',
   '\0': r'\0',
   '\1': r'\1',
   '\2': r'\2',
   '\3': r'\3',
   '\4': r'\4',
   '\5': r'\5',
   '\6': r'\6',
   '\7': r'\7',
   '\8': r'\8',
   '\9': r'\9'
}


def raw(text):
    """Returns a raw string representation of text"""
    new_str = ''
    for char in text:
        try:
            new_str += trans_map[char]
        except KeyError:
            new_str += char

    return new_str


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
        if "smiles2k" in self._get_all_tables():
            print("ERROR: Table smiles2k already exists.")
            return
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
        if "results2k" in self._get_all_tables():
            print("ERROR: Table results2k already exists.")
            return
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
        if "results8m" in self._get_all_tables():
            print("ERROR: Table results8m already exists.")
            return
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
        if "smiles8m" in self._get_all_tables():
            print("ERROR: Table smiles8m already exists.")
            return
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
            print("ERROR: Table {} not found in database.".format(table_name))
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
            print("ERROR: Source file not found in local.")
            return

        if table_name not in self._tables:
            print("ERROR: Table {} not found in database.".format(table_name))
            return

        with open(f, 'r') as ff:
            next(ff)
            self._cursor.copy_from(ff, table_name, sep=',', )
        self._connection.commit()

    def _create_distinct_table(self, table_name='smiles2k'):
        """
        create table with distinct canonical smiles strings
        :param table_name:
        :return:
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database.".format(table_name))
            return

        if "distinct_{}".format(table_name) in self._tables:
            print("ERROR: Table distinct_{} already exists.".format(table_name))
            return

        if "tmp_{}".format(table_name) in self._tables:
            print("ERROR: Table tmp_{} already exists.".format(table_name))
            return

        query = '''
        CREATE TABLE tmp_{}(
        id integer,
        canonical_smiles text,
        molwt float
        );
        '''.format(table_name)
        self._cursor.execute(query)
        self._connection.commit()

        err = []
        count = 0
        count_chunk = 0
        print(count_chunk + 1)
        for i in range(self.get_num_rows(table_name)):
            if count == 1000:
                count_chunk += 1
                print(count_chunk + 1)
                count = 0
                self._connection.commit()
            count += 1
            id_useless, sm, mw = self._fetch_row_by_index(i, table_name)
            try:
                cs = Chem.CanonSmiles(raw(sm))
            except:
                err.append(i)
                continue
            sql = "INSERT INTO tmp_{}".format(table_name)
            sql = sql + " VALUES(%s, %s, %s)"
            self._cursor.execute(sql, (i, cs, mw))
        self._connection.commit()
        print("error list:")
        print(err)
        allq = '''
        SELECT * INTO distinct_{} FROM(
        SELECT MIN(id) AS id,canonical_smiles,MAX(molwt) as molwt  FROM tmp_{} GROUP BY canonical_smiles
         ) tmp;
        '''.format(table_name, table_name)
        self._cursor.execute(allq)
        self._connection.commit()
        self._drop_table('tmp_{}'.format(table_name))

    def _create_err_table(self, table_name):
        """
        create a table for logging run-time-errored canonical smiles
        private cuz only done once
        :param table_name: name of table
        :return:
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database.".format(table_name))
            return

        if "err_{}".format(table_name) in self._tables:
            print("ERROR: Table err_{} already exists.".format(table_name))
            return

        query = '''
        CREATE TABLE err_{}(
        id integer,
        canonical_smiles text,
        molwt float); 
        '''.format(table_name)
        self._cursor.execute(query)
        self._connection.commit()
        self._tables = [i[0] for i in self._get_all_tables()]

    def _create_initial_mol(self, table_name='smiles2k'):
        """
        generate initial molecule table from smiles table
        private cuz only done once for each table
        :param table_name: the name of table to be generated from
        :return:
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database".format(table_name))
            return

        query = '''
        SELECT * INTO mols_{} FROM (SELECT id, mol_from_smiles(smiles::cstring) m FROM {}) tmp WHERE m IS NOT NULL;
        '''.format(table_name, table_name)
        self._cursor.execute(query)
        self._connection.commit()
        self._tables = [i[0] for i in self._get_all_tables()]

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

    def _fetch_row_by_index(self, i, table_name='smiles2k'):
        """
        fetch row from smiles table by index
        :param i: index
        :param table_name: smiles table to be queried
        :return: smiles
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database.".format(table_name))
            return

        query = '''
        SELECT * FROM {} WHERE id = {}
        '''.format(table_name, i)
        self._cursor.execute(query)
        result = self._cursor.fetchone()
        return result

    def _fetch_canonical_smiles_by_index(self, i, table_name='distinct_smiles2k'):
        """
        fetch canonical smiles from distinct smiles table by index
        :param i: index
        :param table_name: distinct smiles table to be queried
        :return: canonical smiles
        """
        query = '''
        SELECT canonical_smiles FROM {} WHERE id = {}
        '''.format(table_name, i)
        self._cursor.execute(query)
        result = self._cursor.fetchone()
        if len(result) == 0:
            return None
        else:
            return raw(result[0])

    def _fetch_smiles_by_index(self, i, table_name='smiles2k'):
        """
        fetch smiles from smiles table by index
        :param i: index
        :param table_name: smiles table to be queried
        :return: smiles
        """
        query = '''
        SELECT smiles FROM {} WHERE id = {}
        '''.format(table_name, i)
        self._cursor.execute(query)
        result = self._cursor.fetchone()
        if len(result) == 0:
            return None
        else:
            return raw(result[0])

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
        result = self._cursor.fetchone()
        return result[0]

    def fetch_ids(self, table_name):
        """
        fetch id list
        :param table_name: table to be queried
        :return:
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database".format(table_name))
            return
        query = '''
        SELECT id FROM {}
        '''.format(table_name)
        self._cursor.execute(query)
        result = [i[0] for i in self._cursor.fetchall()]
        result.sort(key=None, reverse=False)
        return result

    def fetch_canonical_smiles_by_index(self, i, table_name='distinct_smiles2k'):
        """
        fetch canonical smiles from distinct smiles table by index
        :param i: index
        :param table_name: distinct smiles table to be queried
        :return: canonical smiles
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database".format(table_name))
            return

        query = '''
        SELECT canonical_smiles FROM {} WHERE id = {}
        '''.format(table_name, i)
        self._cursor.execute(query)
        result = self._cursor.fetchone()
        if len(result) == 0:
            return None
        else:
            return raw(result[0])

    def fetch_smiles_by_index(self, i, table_name='smiles2k'):
        """
        fetch smiles from smiles table by index
        :param i: index
        :param table_name: smiles table to be queried
        :return: smiles
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database".format(table_name))
            return

        query = '''
        SELECT smiles FROM {} WHERE id = {}
        '''.format(table_name, i)
        self._cursor.execute(query)
        result = self._cursor.fetchone()
        if len(result) == 0:
            return None
        else:
            return raw(result[0])

    def dbtest(self):
        """
        for personal tests
        :return:
        """
        # self.create_err_table("distinct_smiles2k")
        # self.insert_error(89, "distinct_smiles2k")
        self._drop_table("err_distinct_smiles2k")

    def create_err_table(self, table_name):
        """

        :param table_name:
        :return:
        """
        self._create_err_table(table_name)

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

    def insert_error(self, err_id, table_name):
        """
        insert an error smiles record to error table
        :param err_id:
        :param table_name:
        :return:
        """
        if table_name not in self._tables:
            print("ERROR: Table {} not found in database".format(table_name))
            return
        if "err_{}".format(table_name) not in self._tables:
            print("ERROR: Table err_{} not found in database".format(table_name))
            return

        rec = self._fetch_row_by_index(err_id, table_name)
        sql = "INSERT INTO err_{}".format(table_name)
        sql = sql + " VALUES(%s, %s, %s)"
        self._cursor.execute(sql, rec)
        self._connection.commit()

    def reset_results2k(self):
        """

        :return:
        """
        self._drop_table('results2k')
        self._create_table_results2k()
        self._drop_table('err_distinct_smiles2k')

    def reset_results8m(self):
        """

        :return:
        """
        self._drop_table('results8m')
        self._create_table_results8m()
        self._drop_table('err_distinct_smiles8m')

    def get_err_ids(self, table_name):
        """
        get all ids that show invalid prediction
        :param table_name:
        :return:
        """
        if table_name not in self._tables:
            print('ERROR: table name not found in database')
            return

        query = '''
        SELECT id FROM err_{};
        '''.format(table_name)
        self._cursor.execute(query)
        result = self._cursor.fetchall()
        print(result)


if __name__ == '__main__':
    dbu = DBUtil()
    dbu.reset_results2k()
    # dbu.dbtest()
    # dbu.get_err_ids('results2k')
