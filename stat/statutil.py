"""
platform: win
env: any with psycopg and rdkit
name: statutil.py
basic statistics
"""
from dbutil.dbutil import *


class ToyStat:
    _table_name = None
    _db_helper = None
    _min_pk = None
    _max_pk = None

    def _maximum_pk(self):
        """
        get the maximum value of pk
        :return:
        """
        return self._db_helper.get_maximum_pk(self._table_name)

    def _minimum_pk(self):
        """
        get the minimum value of pk
        :return:
        """
        return self._db_helper.get_minimum_pk(self._table_name)

    def __init__(self, table_name="results8m"):
        self._table_name = table_name
        self._db_helper = DBUtil()
        self._min_pk = self._minimum_pk()
        self._max_pk = self._maximum_pk()

    # callables
    def total_num_collected(self):
        """
        get the total number of collected items
        :return:
        """
        return self._db_helper.get_num_rows(self._table_name)

    def num_with_bounded_pk(self, lower_bound=1.0, upper_bound=11.0):
        """
        get the number of collected items with two pk boundaries
        :param lower_bound: lower pk boundary
        :param upper_bound: upper pk boundary
        :return:
        """
        if lower_bound > upper_bound:
            print("Error: lower bound greater than upper bound")
            return 0
        if lower_bound > self._max_pk:
            return 0
        if upper_bound < self._min_pk:
            return 0
        return self._db_helper.count_aff_with_bounded_pk(self._table_name, lower_bound, upper_bound)

    def num_greater_than(self, cutoff=1.0):
        """
        get the number of collected items greater than a cutoff
        :param cutoff:
        :return:
        """
        if cutoff > self._max_pk:
            return 0
        return self._db_helper.count_aff_greater_than(self._table_name, cutoff)

    def num_lesser_than(self, cutoff=11.0):
        """
        get the number of collected items lesser than a cutoff
        :param cutoff:
        :return:
        """
        if cutoff < self._min_pk:
            return 0
        return self._db_helper.count_aff_lesser_than(self._table_name, cutoff)

    def maximum_pk(self):
        """
        get maximum pk value
        :return:
        """
        return self._max_pk

    def minimum_pk(self):
        """
        get minimum pk value
        :return:
        """
        return self._min_pk

    def toy_stat_test(self):
        """
        class tester
        :return:
        """
        pass


def test():
    stat_helper = ToyStat()
    print(stat_helper.num_greater_than(10.2))


if __name__ == '__main__':
    test()
