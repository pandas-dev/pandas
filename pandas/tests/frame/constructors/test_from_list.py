
import pytest
import pandas as pd
import pandas._testing as tm
import numpy as np



class TestFromList:
     
    #  passing lists to DataFrame.__init__
    # related GH#6140
    def test_start_with_a_List_and_directions_parameter_is_not_setted(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        dict = {'a': 1, 'b': 2,'c':3,'d':4}
        result = pd.DataFrame(data=list_data,columns=list_columns,index=[0])
        expected  = pd.DataFrame(dict,index=[0])
        tm.assert_frame_equal(result,expected)

    def test_start_with_a_List_and_directions_parameter_is_not_setted_and_index_parameter_is_not_setted(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        dict = {'a': 1, 'b': 2,'c':3,'d':4}
        result = pd.DataFrame(data=list_data,columns=list_columns)
        expected  = pd.DataFrame(dict,index=[0])
        tm.assert_frame_equal(result,expected)
    
    
    def test_direction_parameter_len_higher_then_two_and_direction_setted_as_arrow(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[0],direction=["arrow",None,"aleatory_variable"])
           
    def test_direction_parameter_len_lower_then_two_and_direction_setted_as_arrow(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[0],direction=["arrow"])
            
    
    def test_direction_parameter_len_higher_then_two_and_direction_setted_as_columns(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[0],direction=["columns",None,"aleatory_variable"])

    def test_direction_parameter_len_lower_then_two_and_direction_setted_as_columns(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[0],direction=["columns"])

    
    def test_direction_parameter_setted_as_arrow_with_second_element_different_from_None(self): 
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[0],direction=["arrow","Aleatory_variable"])

    def test_create_a_dataframe_when_direction_parameter_is_setted_as_columns(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        nan_list = list()
        for i in range(0,4):
            nan_list.append(np.nan)


        dict = {'a': [1,2,3,4], 'b':nan_list,'c':nan_list,'d':nan_list}
        result = pd.DataFrame(data=list_data,columns=list_columns,direction=["columns","a"])
        expected = pd.DataFrame(data=dict, index=[0, 1, 2, 3])
        tm.assert_frame_equal(result,expected)


    def test_direction_parameter_columns_with_second_element_None(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[0],direction=["columns", None])


    def test_direction_parameter_setted_as_columns_and_index_parameter_is_setted(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[2],direction=["columns","a"])

    def test_invalid_values_in_direction_parameter(self):
        list_data = [1,2,3,4]
        list_columns = ["a","b","c","d"]
        with pytest.raises(Exception):
            pd.DataFrame(data=list_data,columns=list_columns,index=[2],direction=["Invalid",5])
    

    