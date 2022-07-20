"""
platform: win
env: any but with psycopg2 lib
name: dbconn.py
database: you have to create it first locally
enable connection to db and some basic operations
"""
import os.path
from config_win import *
import psycopg2


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
        query = '''create table smiles2k(
        id serial primary key,
        smiles text,
        molwt float); 
        '''
        self._cursor.execute(query)
        self._connection.commit()

    def _create_table_smiles8m(self):
        """
        create a table containing 8m smiles strings
        for test
        private cuz only done once
        :return:
        """
        query = '''
        create table smiles8m(
        id serial primary key,
        smiles text,
        molwt float); 
        '''
        self._cursor.execute(query)
        self._connection.commit()



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
        SELECT * INTO mols_{} FROM (select id, mol_from_smiles(smiles::cstring) m from {}) tmp WHERE m IS NOT NULL;
        '''.format(table_name, table_name)
        self._cursor.execute(query)
        self._connection.commit()

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

    def dbtest(self):
        """
        for personal tests
        :return:
        """
        self._create_initial_mol('smiles2k')


if __name__ == '__main__':
    dbu = DBUtil()
    dbu.dbtest()
