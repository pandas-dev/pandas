"""
Test cases for panda to sql
"""
# pylint: disable=broad-except
import os
from pathlib import Path
import numpy as np
from datetime import datetime, date
from pandas.core.frame import DataFrame
from pandas.core.reshape.concat import concat
from pandas.core.reshape.merge import merge
from pandas.io.parsers import read_csv
from freezegun import freeze_time
from pandas.core.sql.sql_to_pandas import SqlToPandas
from pandas.core.sql.sql_exception import InvalidQueryException, DataFrameDoesNotExist
from pandas.core.sql.sql_to_pandas import register_temp_table
from pandas.util.testing import assert_frame_equal


DATA_PATH = os.path.join(Path(__file__).parent, "data")


# Import the data for testing
FOREST_FIRES = read_csv(os.path.join(DATA_PATH, 'forestfires.csv'))
DIGIMON_MON_LIST = read_csv(os.path.join(DATA_PATH, 'DigiDB_digimonlist.csv'))
DIGIMON_MOVE_LIST = read_csv(os.path.join(DATA_PATH, 'DigiDB_movelist.csv'))
DIGIMON_SUPPORT_LIST = read_csv(os.path.join(DATA_PATH, 'DigiDB_supportlist.csv'))

# Name change is for name interference
DIGIMON_MON_LIST['mon_attribute'] = DIGIMON_MON_LIST['Attribute']
DIGIMON_MOVE_LIST['move_attribute'] = DIGIMON_MOVE_LIST['Attribute']


def register_env_tables():
    """
    Returns all globals but in lower case
    :return:
    """
    for variable_name in globals():
        variable = globals()[variable_name]
        if isinstance(variable, DataFrame):
            register_temp_table(frame=variable, table_name=variable_name)

register_env_tables()


def sql_to_pandas_with_vars(sql: str):
    """
    Preset with data in SqlToPandas class
    :param sql: Sql query
    :return: SqlToPandasClass with
    """
    return SqlToPandas(sql).data_frame


def test_for_valid_query():
    """
    Test that exception is raised for invalid query
    :return:
    """
    sql = "hello world!"
    try:
        sql_to_pandas_with_vars(sql)
    except InvalidQueryException as err:
        assert isinstance(err, InvalidQueryException)


def test_select_star():
    """
    Tests the simple select * case
    :return:
    """
    myframe = sql_to_pandas_with_vars("select * from forest_fires")
    assert FOREST_FIRES.equals(myframe)


def test_case_insensitivity():
    """
    Tests to ensure that the sql is case insensitive for table names
    :return:
    """
    assert FOREST_FIRES.equals(sql_to_pandas_with_vars("select * from FOREST_fires"))


def test_select_specific_fields():
    """
    Tests selecting specific fields
    :return:
    """
    myframe = sql_to_pandas_with_vars("select temp, RH, wind, rain as water, area from forest_fires")
    pandas_frame = FOREST_FIRES[['temp', 'RH', 'wind', 'rain', 'area']].rename(columns={'rain': 'water'})
    assert myframe.equals(pandas_frame)


def test_type_conversion():
    """
    Tests sql as statements
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select cast(temp as int64), cast(RH as float64) my_rh, wind, rain, area,
    cast(2.0 as int64) my_int, cast(3 as float64) as my_float, cast(7 as object) as my_object, 
    cast(0 as bool) as my_bool from forest_fires""")
    fire_frame = FOREST_FIRES[['temp', 'RH', 'wind', 'rain', 'area']].rename(columns={'RH': 'my_rh'})
    fire_frame["my_int"] = 2
    fire_frame["my_float"] = 3
    fire_frame["my_object"] = str(7)
    fire_frame["my_bool"] = 0
    pandas_frame = fire_frame.astype({'temp': 'int64', 'my_rh': 'float64', 'my_int': 'int64', 'my_float': 'float64',
                                      'my_bool': 'bool'})
    assert pandas_frame.equals(my_frame)


def test_for_non_existent_table():
    """
    Check that exception is raised if table does not exist
    :return:
    """
    try:
        sql_to_pandas_with_vars("select * from a_table_that_is_not_here")
    except Exception as err:
        assert isinstance(err, DataFrameDoesNotExist)


def test_using_math():
    """
    Test the mathematical operations and order of operations
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select temp, 1 + 2 * 3 as my_number from forest_fires")
    pandas_frame = FOREST_FIRES[['temp']].copy()
    pandas_frame['my_number'] = 1 + 2 * 3
    assert pandas_frame.equals(my_frame)


def test_distinct():
    """
    Test use of the distinct keyword
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select distinct area, rain from forest_fires")
    pandas_frame = FOREST_FIRES[['area', 'rain']].copy()
    pandas_frame.drop_duplicates(keep='first', inplace=True)
    pandas_frame.reset_index(inplace=True)
    pandas_frame.drop(columns='index', inplace=True)
    assert pandas_frame.equals(my_frame)


def test_subquery():
    """
    Test ability to perform subqueries
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select * from (select area, rain from forest_fires) rain_area")
    pandas_frame = FOREST_FIRES[['area', 'rain']].copy()
    assert pandas_frame.equals(my_frame)


def test_join_no_inner():
    """
    Test join
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list join
            digimon_move_list
            on digimon_mon_list.attribute = digimon_move_list.attribute""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, on="Attribute")
    assert merged_frame.equals(my_frame)


def test_join_wo_specifying_table():
    """
    Test join where table isn't specified in join
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """
        select * from digimon_mon_list join
        digimon_move_list
        on mon_attribute = move_attribute
        """)
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, left_on="mon_attribute", right_on="move_attribute")
    assert merged_frame.equals(my_frame)


def test_join_w_inner():
    """
    Test join
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list inner join
            digimon_move_list
            on digimon_mon_list.attribute = digimon_move_list.attribute""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, on="Attribute")
    assert merged_frame.equals(my_frame)


def test_outer_join_no_outer():
    """
    Test outer join
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list full outer join
            digimon_move_list
            on digimon_mon_list.type = digimon_move_list.type""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, how="outer", on="Type")
    assert merged_frame.equals(my_frame)


def test_outer_join_w_outer():
    """
    Test outer join
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list full join
            digimon_move_list
            on digimon_mon_list.type = digimon_move_list.type""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, how="outer", on="Type")
    assert merged_frame.equals(my_frame)


def test_left_joins():
    """
    Test right, left, inner, and outer joins
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list left join
            digimon_move_list
            on digimon_mon_list.type = digimon_move_list.type""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, how="left", on="Type")
    assert merged_frame.equals(my_frame)


def test_left_outer_joins():
    """
    Test right, left, inner, and outer joins
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list left outer join
            digimon_move_list
            on digimon_mon_list.type = digimon_move_list.type""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, how="left", on="Type")
    assert merged_frame.equals(my_frame)


def test_right_joins():
    """
    Test right, left, inner, and outer joins
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list right join
            digimon_move_list
            on digimon_mon_list.type = digimon_move_list.type""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, how="right", on="Type")
    assert merged_frame.equals(my_frame)


def test_right_outer_joins():
    """
    Test right, left, inner, and outer joins
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list right outer join
            digimon_move_list
            on digimon_mon_list.type = digimon_move_list.type""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, how="right", on="Type")
    assert merged_frame.equals(my_frame)


def test_cross_joins():
    """
    Test right, left, inner, and outer joins
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from digimon_mon_list cross join
            digimon_move_list
            on digimon_mon_list.type = digimon_move_list.type""")
    pandas_frame1 = DIGIMON_MON_LIST
    pandas_frame2 = DIGIMON_MOVE_LIST
    merged_frame = pandas_frame1.merge(pandas_frame2, how="outer", on="Type")
    assert merged_frame.equals(my_frame)


def test_group_by():
    """
    Test group by constraint
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select month, day from forest_fires group by month, day""")
    pandas_frame = FOREST_FIRES.groupby(["month", "day"]).size().to_frame('size').reset_index().drop(columns=['size'])
    assert pandas_frame.equals(my_frame)


def test_avg():
    """
    Test the avg
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select avg(temp) from forest_fires")
    pandas_frame = FOREST_FIRES.agg({'temp': np.mean}).to_frame('mean_temp').reset_index().drop(columns=['index'])
    assert pandas_frame.equals(my_frame)

def test_sum():
    """
    Test the sum
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select sum(temp) from forest_fires")
    pandas_frame = FOREST_FIRES.agg({'temp': np.sum}).to_frame('sum_temp').reset_index().drop(columns=['index'])
    assert pandas_frame.equals(my_frame)


def test_max():
    """
    Test the max
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select max(temp) from forest_fires")
    pandas_frame = FOREST_FIRES.agg({'temp': np.max}).to_frame('max_temp').reset_index().drop(columns=['index'])
    assert pandas_frame.equals(my_frame)

def test_min():
    """
    Test the min
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select min(temp) from forest_fires")
    pandas_frame = FOREST_FIRES.agg({'temp': np.min}).to_frame('min_temp').reset_index().drop(columns=['index'])
    assert pandas_frame.equals(my_frame)

def test_multiple_aggs():
    """
    Test multiple aggregations
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select min(temp), max(temp), avg(temp), max(wind) from forest_fires")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame['min_temp'] = FOREST_FIRES.temp.copy()
    pandas_frame['max_temp'] = FOREST_FIRES.temp.copy()
    pandas_frame['mean_temp'] = FOREST_FIRES.temp.copy()
    pandas_frame = pandas_frame.agg({'min_temp': np.min, 'max_temp': np.max, 'mean_temp': np.mean, 'wind': np.max})
    pandas_frame.rename({'wind': 'max_wind'}, inplace=True)
    pandas_frame = pandas_frame.to_frame().transpose()
    assert pandas_frame.equals(my_frame)


def test_agg_w_groupby():
    """
    Test using aggregates and group by together
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select day, month, min(temp), max(temp) from forest_fires group by day, month")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame['min_temp'] = pandas_frame.temp
    pandas_frame['max_temp'] = pandas_frame.temp
    pandas_frame = pandas_frame.groupby(["day", "month"]).aggregate({'min_temp': np.min, 'max_temp': np.max})\
        .reset_index()
    assert pandas_frame.equals(my_frame)


def test_where_clause():
    """
    Test where clause
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select * from forest_fires where month = 'mar'""")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame = pandas_frame[pandas_frame.month == 'mar'].reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_order_by():
    """
    Test order by clause
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select * from forest_fires order by temp desc, wind asc, area""")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame.sort_values(by=['temp', 'wind', 'area'], ascending=[0, 1, 1], inplace=True)
    pandas_frame.reset_index(drop=True, inplace=True)
    assert pandas_frame.equals(my_frame)


def test_limit():
    """
    Test limit clause
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select * from forest_fires limit 10""")
    pandas_frame = FOREST_FIRES.copy().head(10)
    assert pandas_frame.equals(my_frame)


def test_having():
    """
    Test having clause
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select min(temp) from forest_fires having min(temp) > 2")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame['min_temp'] = FOREST_FIRES['temp']
    aggregated_df = pandas_frame.aggregate({'min_temp': 'min'}).to_frame().transpose()
    pandas_frame = aggregated_df[aggregated_df['min_temp'] > 2]
    assert pandas_frame.equals(my_frame)


def test_having_with_group_by():
    """
    Test having clause
    :return:
    """
    my_frame = sql_to_pandas_with_vars("select day, min(temp) from forest_fires group by day having min(temp) > 5")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame['min_temp'] = FOREST_FIRES['temp']
    pandas_frame = pandas_frame[['day', 'min_temp']].groupby('day').aggregate({'min_temp': np.min})
    pandas_frame = pandas_frame[pandas_frame['min_temp'] > 5].reset_index()
    assert pandas_frame.equals(my_frame)


def test_operations_between_columns_and_numbers():
    """
    Tests operations between columns
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select temp * wind + rain / dmc + 37 from forest_fires""")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame['temp_mul_wind_add_rain_div_dmc_add_37'] = pandas_frame['temp'] * pandas_frame['wind'] + \
                                                            pandas_frame['rain'] / pandas_frame['DMC'] + 37
    pandas_frame = pandas_frame['temp_mul_wind_add_rain_div_dmc_add_37'].to_frame()
    assert pandas_frame.equals(my_frame)


def test_select_star_from_multiple_tables():
    """
    Test selecting from two different tables
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select * from forest_fires, digimon_mon_list""")
    forest_fires = FOREST_FIRES.copy()
    digimon_mon_list_new = DIGIMON_MON_LIST.copy()
    forest_fires['_temp_id'] = 1
    digimon_mon_list_new['_temp_id'] = 1
    pandas_frame = merge(forest_fires, digimon_mon_list_new, on='_temp_id').drop(columns=['_temp_id'])
    assert pandas_frame.equals(my_frame)


def test_select_columns_from_two_tables_with_same_column_name():
    """
    Test selecting tables
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""select * from forest_fires table1, forest_fires table2""")
    table1 = FOREST_FIRES.copy()
    table2 = FOREST_FIRES.copy()
    table1['_temp_id'] = 1
    table2['_temp_id'] = 1
    pandas_frame = merge(table1, table2, on='_temp_id').drop(columns=['_temp_id'])
    assert pandas_frame.equals(my_frame)


def test_maintain_case_in_query():
    """
    Test nested subqueries
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select wind, rh from forest_fires""")
    pandas_frame = FOREST_FIRES.copy()[['wind', 'RH']].rename(columns={"RH": "rh"})
    assert pandas_frame.equals(my_frame)


def test_nested_subquery():
    """
    Test nested subqueries
    :return:
    """
    my_frame = sql_to_pandas_with_vars(
        """select * from
            (select wind, rh from
              (select * from forest_fires) fires) wind_rh""")
    pandas_frame = FOREST_FIRES.copy()[['wind', 'RH']].rename(columns={"RH": "rh"})
    assert pandas_frame.equals(my_frame)


def test_union():
    """
    Test union in queries
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select * from forest_fires order by wind desc limit 5
     union 
    select * from forest_fires order by wind asc limit 5
    """)
    pandas_frame1 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(5)
    pandas_frame2 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[True]).head(5)
    pandas_frame = concat([pandas_frame1, pandas_frame2], ignore_index=True).drop_duplicates().reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_union_distinct():
    """
    Test union distinct in queries
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
        select * from forest_fires order by wind desc limit 5
         union distinct
        select * from forest_fires order by wind asc limit 5
        """)
    pandas_frame1 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(5)
    pandas_frame2 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[True]).head(5)
    pandas_frame = concat([pandas_frame1, pandas_frame2], ignore_index=True).drop_duplicates().reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_union_all():
    """
    Test union distinct in queries
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
        select * from forest_fires order by wind desc limit 5
         union all
        select * from forest_fires order by wind asc limit 5
        """)
    pandas_frame1 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(5)
    pandas_frame2 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[True]).head(5)
    pandas_frame = concat([pandas_frame1, pandas_frame2], ignore_index=True).reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_intersect_distinct():
    """
    Test union distinct in queries
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
            select * from forest_fires order by wind desc limit 5
             intersect distinct
            select * from forest_fires order by wind desc limit 3
            """)
    pandas_frame1 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(5)
    pandas_frame2 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(3)
    pandas_frame = merge(left=pandas_frame1, right=pandas_frame2, how='inner', on=list(pandas_frame1.columns))
    assert pandas_frame.equals(my_frame)


def test_except_distinct():
    """
    Test except distinct in queries
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
                select * from forest_fires order by wind desc limit 5
                 except distinct
                select * from forest_fires order by wind desc limit 3
                """)
    pandas_frame1 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(5)
    pandas_frame2 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(3)
    pandas_frame = pandas_frame1[~pandas_frame1.isin(pandas_frame2).all(axis=1)].drop_duplicates().reset_index(
        drop=True)
    assert pandas_frame.equals(my_frame)


def test_except_all():
    """
    Test except distinct in queries
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
                select * from forest_fires order by wind desc limit 5
                 except all
                select * from forest_fires order by wind desc limit 3
                """)
    pandas_frame1 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(5)
    pandas_frame2 = FOREST_FIRES.copy().sort_values(by=['wind'], ascending=[False]).head(3)
    pandas_frame = pandas_frame1[~pandas_frame1.isin(pandas_frame2).all(axis=1)].reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_between_operator():
    """
    Test using between operator
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select * from forest_fires
    where wind between 5 and 6
    """)
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame = pandas_frame[(pandas_frame.wind >= 5) & (pandas_frame.wind <= 6)].reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_in_operator():
    """
    Test using in operator in a sql query
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select * from forest_fires where day in ('fri', 'sun')
    """)
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame = pandas_frame[pandas_frame.day.isin(('fri', 'sun'))].reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_in_operator_expression_numerical():
    """
    Test using in operator in a sql query
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select * from forest_fires where X + 1 in (5, 9)
    """)
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame = pandas_frame[(pandas_frame['X'] + 1).isin((5, 9))].reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_not_in_operator():
    """
    Test using in operator in a sql query
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select * from forest_fires where day not in ('fri', 'sun')
    """)
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame = pandas_frame[~pandas_frame.day.isin(('fri', 'sun'))].reset_index(drop=True)
    assert pandas_frame.equals(my_frame)


def test_case_statement_w_name():
    """
    Test using case statements
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
        select case when wind > 5 then 'strong' when wind = 5 then 'mid' else 'weak' end as wind_strength from forest_fires
        """)
    pandas_frame = FOREST_FIRES.copy()[['wind']]
    pandas_frame.loc[pandas_frame.wind > 5, 'wind_strength'] = 'strong'
    pandas_frame.loc[pandas_frame.wind == 5, 'wind_strength'] = 'mid'
    pandas_frame.loc[~((pandas_frame.wind == 5) | (pandas_frame.wind > 5)), 'wind_strength'] = 'weak'
    pandas_frame.drop(columns=['wind'], inplace=True)
    assert pandas_frame.equals(my_frame)


def test_case_statement_w_no_name():
    """
    Test using case statements
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
        select case when wind > 5 then 'strong' when wind = 5 then 'mid' else 'weak' end from forest_fires
        """)
    pandas_frame = FOREST_FIRES.copy()[['wind']]
    pandas_frame.loc[pandas_frame.wind > 5, '_expression0'] = 'strong'
    pandas_frame.loc[pandas_frame.wind == 5, '_expression0'] = 'mid'
    pandas_frame.loc[~((pandas_frame.wind == 5) | (pandas_frame.wind > 5)), '_expression0'] = 'weak'
    pandas_frame.drop(columns=['wind'], inplace=True)
    assert pandas_frame.equals(my_frame)


def test_case_statement_w_other_columns_as_reult():
    """
    Test using case statements
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
        select case when wind > 5 then month when wind = 5 then 'mid' else day end from forest_fires
        """)
    pandas_frame = FOREST_FIRES.copy()[['wind']]
    pandas_frame.loc[pandas_frame.wind > 5, '_expression0'] = FOREST_FIRES['month']
    pandas_frame.loc[pandas_frame.wind == 5, '_expression0'] = 'mid'
    pandas_frame.loc[~((pandas_frame.wind == 5) | (pandas_frame.wind > 5)), '_expression0'] = FOREST_FIRES['day']
    pandas_frame.drop(columns=['wind'], inplace=True)
    assert pandas_frame.equals(my_frame)


def test_rank_statement_one_column():
    """
    Test rank statement
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select wind, rank() over(order by wind) as wind_rank from forest_fires
    """)
    pandas_frame = FOREST_FIRES.copy()[['wind']]
    pandas_frame['wind_rank'] = pandas_frame.wind.rank(method='min').astype('int')
    assert pandas_frame.equals(my_frame)


def test_rank_statement_many_columns():
    """
    Test rank statement
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select wind, rain, month, rank() over(order by wind desc, rain asc, month) as rank from forest_fires
    """)
    pandas_frame = FOREST_FIRES.copy()[['wind', 'rain', 'month']]
    pandas_frame.sort_values(by=['wind', 'rain', 'month'], ascending=[False, True, True], inplace=True)
    pandas_frame.reset_index(inplace=True)
    rank_map = {}
    rank_counter = 1
    rank_offset = 0
    pandas_frame['rank'] = 0
    rank_series = pandas_frame['rank'].copy()
    for row_num, row in enumerate(pandas_frame.iterrows()):
        key = "".join(map(str, list(list(row)[1])[1:4]))
        if rank_map.get(key):
            rank_offset += 1
            rank = rank_map[key]
        else:
            rank = rank_counter + rank_offset
            rank_map[key] = rank
            rank_counter += 1
        rank_series[row_num] = rank
    pandas_frame['rank'] = rank_series
    pandas_frame.set_index('index', inplace=True)
    pandas_frame.sort_index(inplace=True)
    assert pandas_frame.equals(my_frame)


def test_dense_rank_statement_many_columns():
    """
    Test dense_rank statement
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select wind, rain, month, dense_rank() over(order by wind desc, rain asc, month) as rank from forest_fires
    """)
    pandas_frame = FOREST_FIRES.copy()[['wind', 'rain', 'month']]
    pandas_frame.sort_values(by=['wind', 'rain', 'month'], ascending=[False, True, True], inplace=True)
    pandas_frame.reset_index(inplace=True)
    rank_map = {}
    rank_counter = 1
    pandas_frame['rank'] = 0
    rank_series = pandas_frame['rank'].copy()
    for row_num, row in enumerate(pandas_frame.iterrows()):
        key = "".join(map(str, list(list(row)[1])[1:4]))
        if rank_map.get(key):
            rank = rank_map[key]
        else:
            rank = rank_counter
            rank_map[key] = rank
            rank_counter += 1
        rank_series[row_num] = rank
    pandas_frame['rank'] = rank_series
    pandas_frame.set_index('index', inplace=True)
    pandas_frame.sort_index(inplace=True)
    assert pandas_frame.equals(my_frame)


def test_rank_over_partition_by():
    """
    Test rank partition by statement
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select wind, rain, month, day, rank() over(partition by day order by wind desc, rain asc, month) as rank
    from forest_fires
    """)
    pandas_frame = FOREST_FIRES.copy()[['wind', 'rain', 'month', 'day']]
    partition_slice = 4
    rank_map = {}
    partition_rank_counter = {}
    partition_rank_offset = {}
    pandas_frame.sort_values(by=['wind', 'rain', 'month'], ascending=[False, True, True], inplace=True)
    pandas_frame.reset_index(inplace=True)
    pandas_frame['rank'] = 0
    rank_series = pandas_frame['rank'].copy()
    for row_num, series_tuple in enumerate(pandas_frame.iterrows()):
        row = series_tuple[1]
        row_list = list(row)[1:partition_slice]
        partition_list = list(row)[partition_slice:5]
        key = str(row_list)
        partition_key = str(partition_list)
        if rank_map.get(partition_key):
            if rank_map[partition_key].get(key):
                partition_rank_counter[partition_key] += 1
                rank = rank_map[partition_key][key]
            else:
                partition_rank_counter[partition_key] += 1
                rank = partition_rank_counter[partition_key] + partition_rank_offset[partition_key]
                rank_map[partition_key][key] = rank
        else:
            rank = 1
            rank_map[partition_key] = {}
            partition_rank_counter[partition_key] = 1
            partition_rank_offset[partition_key] = 0
            rank_map[partition_key][key] = rank
        rank_series[row_num] = rank
    pandas_frame['rank'] = rank_series
    pandas_frame.set_index('index', inplace=True)
    pandas_frame.sort_index(inplace=True)
    assert pandas_frame.equals(my_frame)

def test_dense_rank_over_partition_by():
    """
    Test rank partition by statement
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select wind, rain, month, day, dense_rank() over(partition by day order by wind desc, rain asc, month) as rank
    from forest_fires
    """)
    pandas_frame = FOREST_FIRES.copy()[['wind', 'rain', 'month', 'day']]
    partition_slice = 4
    rank_map = {}
    partition_rank_counter = {}
    pandas_frame.sort_values(by=['wind', 'rain', 'month'], ascending=[False, True, True], inplace=True)
    pandas_frame.reset_index(inplace=True)
    pandas_frame['rank'] = 0
    rank_series = pandas_frame['rank'].copy()
    for row_num, series_tuple in enumerate(pandas_frame.iterrows()):
        row = series_tuple[1]
        row_list = list(row)[1:partition_slice]
        partition_list = list(row)[partition_slice:]
        key = str(row_list)
        partition_key = str(partition_list)
        if rank_map.get(partition_key):
            if rank_map[partition_key].get(key):
                rank = rank_map[partition_key][key]
            else:
                partition_rank_counter[partition_key] += 1
                rank = partition_rank_counter[partition_key]
                rank_map[partition_key][key] = rank
        else:
            rank = 1
            rank_map[partition_key] = {}
            partition_rank_counter[partition_key] = 1
            rank_map[partition_key][key] = rank
        rank_series[row_num] = rank
    pandas_frame['rank'] = rank_series
    pandas_frame.set_index('index', inplace=True)
    pandas_frame.sort_index(inplace=True)
    assert pandas_frame.equals(my_frame)


def test_set_string_value_as_column_value():
    """
    Select a string like 'Yes' as a column value
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select wind, 'yes' as wind_yes from forest_fires""")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame['wind_yes'] = 'yes'
    pandas_frame = pandas_frame[['wind', 'wind_yes']]
    assert pandas_frame.equals(my_frame)


def test_date_cast():
    """
    Select casting a string as a date
    :return:
    """
    my_frame = sql_to_pandas_with_vars("""
    select wind, cast('2019-01-01' as datetime64) as my_date from forest_fires""")
    pandas_frame = FOREST_FIRES.copy()
    pandas_frame['my_date'] = datetime.strptime('2019-01-01', "%Y-%m-%d")
    pandas_frame = pandas_frame[['wind', 'my_date']]
    assert pandas_frame.equals(my_frame)


def test_timestamps():
    """
    Select now() as date
    :return:
    """
    with freeze_time(datetime.now()):
        my_frame = sql_to_pandas_with_vars("""
        select wind, now(), today(), timestamp('2019-01-31', '23:20:32') from forest_fires""")
        pandas_frame = FOREST_FIRES.copy()[['wind']]
        pandas_frame['now()'] = datetime.now()
        pandas_frame['today()'] = date.today()
        pandas_frame['_literal0'] = datetime(2019, 1, 31, 23, 20, 32)
        assert pandas_frame.equals(my_frame)

if __name__ == "__main__":
    pass
