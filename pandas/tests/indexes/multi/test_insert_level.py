import pytest
import pandas as pd
import numpy as np
import pandas._testing as tm



class TestMultiIndexInsertLevel:
    """测试MultiIndex.insert_level方法"""

    def setup_method(self):
        """测试前置准备"""
        # 创建基础测试数据
        self.simple_idx = pd.MultiIndex.from_tuples(
            [('A', 1), ('B', 2), ('C', 3)], names=['level1', 'level2']
        )
        self.empty_idx = pd.MultiIndex.from_tuples([], names=['level1', 'level2'])

    def test_insert_level_basic(self):
        """测试基本功能"""
        # 在位置0插入新层级
        result = self.simple_idx.insert_level(0, 'new_value')
        expected = pd.MultiIndex.from_tuples(
            [('new_value', 'A', 1), ('new_value', 'B', 2), ('new_value', 'C', 3)],
            names=[None, 'level1', 'level2']  # 新插入的层级名称为None
        )
        tm.assert_index_equal(result, expected)

        # 在位置1插入新层级
        result = self.simple_idx.insert_level(1, 'middle')
        expected = pd.MultiIndex.from_tuples(
            [('A', 'middle', 1), ('B', 'middle', 2), ('C', 'middle', 3)],
            names=['level1', None, 'level2']  # 新插入的层级名称为None
        )
        tm.assert_index_equal(result, expected)

    def test_insert_level_with_different_values(self):
        """测试插入不同值的层级"""
        new_values = ['X', 'Y', 'Z']
        result = self.simple_idx.insert_level(1, new_values)
        expected = pd.MultiIndex.from_tuples(
            [('A', 'X', 1), ('B', 'Y', 2), ('C', 'Z', 3)],
            names=['level1', None, 'level2']  # 新插入的层级名称为None
        )
        tm.assert_index_equal(result, expected)

    def test_insert_level_with_name(self):
        """测试指定层级名称"""
        result = self.simple_idx.insert_level(0, 'new_val', name='new_level')
        assert result.names[0] == 'new_level'

    def test_insert_level_edge_positions(self):
        """测试边界位置插入"""
        # 在开始位置插入
        result_start = self.simple_idx.insert_level(0, 'start')
        assert result_start.nlevels == 3

        # 在结束位置插入
        result_end = self.simple_idx.insert_level(2, 'end')
        assert result_end.nlevels == 3

    def test_insert_level_error_cases(self):
        """测试错误情况"""
        # 位置超出范围
        with pytest.raises(ValueError, match="position must be between"):
            self.simple_idx.insert_level(5, 'invalid')

        # 位置为负数
        with pytest.raises(ValueError, match="position must be between"):
            self.simple_idx.insert_level(-1, 'invalid')

        # 值长度不匹配
        with pytest.raises(ValueError, match="Length of values must match"):
            self.simple_idx.insert_level(1, ['too', 'few'])

    def test_insert_level_with_different_data_types(self):
        """测试不同数据类型"""
        # 整数
        result_int = self.simple_idx.insert_level(1, 100)

        # 浮点数
        result_float = self.simple_idx.insert_level(1, 1.5)

        # None值
        result_none = self.simple_idx.insert_level(1, None)

        # 确保都能正常创建
        assert result_int.nlevels == 3
        assert result_float.nlevels == 3
        assert result_none.nlevels == 3

    def test_insert_level_preserves_original(self):
        """测试原索引不被修改"""
        original = self.simple_idx.copy()
        result = self.simple_idx.insert_level(1, 'temp')

        # 原索引应保持不变
        tm.assert_index_equal(original, self.simple_idx)
        # 新索引应有更多层级
        assert result.nlevels == original.nlevels + 1

        def test_debug_names():
            """调试层级名称问题"""
            idx = pd.MultiIndex.from_tuples(
                [('A', 1), ('B', 2), ('C', 3)],
                names=['level1', 'level2']
            )
            print("Original names:", idx.names)

            result = idx.insert_level(0, 'new_value')
            print("Result names:", result.names)

            # 手动创建期望的结果
            expected = pd.MultiIndex.from_tuples(
                [('new_value', 'A', 1), ('new_value', 'B', 2), ('new_value', 'C', 3)],
                names=[None, 'level1', 'level2']  # 注意：新插入的层级名称应该是None
            )
            print("Expected names:", expected.names)