from pandas.arrays import IntervalArray
from pandas import Interval
import pandas._testing as tm

# Unit tests for issue 43629
class TestIntervalArithmetic:
    def test_addition(self):
        cases = [
            # [1, 5] + [2, 7] = [1, 7]
            (Interval(1, 5, 'both'), Interval(2, 7, 'both'), Interval(1, 7, 'both')),
            # (1, 5) + (2, 7) = (1, 7)
            (Interval(1, 5, 'neither'), Interval(2, 7, 'neither'), Interval(1, 7, 'neither')),
            # (1, 5] + (2, 7] = (1, 7]
            (Interval(1, 5, 'right'), Interval(2, 7, 'right'), Interval(1, 7, 'right')),
            # (1, 5] + (5, 7] = (1, 7]
            (Interval(1, 5, 'right'), Interval(5, 7, 'right'), Interval(1, 7, 'right')),
            # (1, 5) + [5, 7] = (1, 7]
            (Interval(1, 5, 'neither'), Interval(5, 7, 'both'), Interval(1, 7, 'right')),
            # (1, 5) + (5, 7) = [(1, 5), (5, 7)]
            (Interval(1, 5, 'neither'), Interval(5, 7, 'neither'), [Interval(1, 5, 'neither'), Interval(5, 7, 'neither')]),
            # (1, 4) + (5, 7] = [(1, 4), (5, 7]]
            (Interval(1, 4, 'neither'), Interval(5, 7, 'right'), [Interval(1, 4, 'neither'), Interval(5, 7, 'right')]),
            # (1, 7) + (5, 6] = (1, 7)
            (Interval(1, 7, 'neither'), Interval(5, 6, 'right'), Interval(1, 7, 'neither')),
            #(1,7] + [1,7) = [1,7]
            (Interval(1, 7, 'right'), Interval(1, 7, 'left'), Interval(1, 7, 'both'))
        ]
        for (interval1, interval2, expected) in cases:
            actual1 = interval1 + interval2
            actual2 = interval2 + interval1
            assert actual1 == expected
            assert actual2 == expected


    def test_multiplication(self):
        cases = [
            # (1, 5) * (2, 7) = (2, 5)
            (Interval(1, 5, 'neither'), Interval(2, 7, 'neither'), Interval(2, 5, 'neither')),
            # (1, 5) * (6, 7) = (0, 0)
            (Interval(1, 5, 'neither'), Interval(6, 7, 'neither'), Interval(0, 0, 'neither')),
            # (1, 5] * (2, 7) = (2, 5]
            (Interval(1, 5, 'right'), Interval(2, 7, 'neither'), Interval(2, 5, 'right')),
            # (1, 7) * [3, 5] = [3, 5]
            (Interval(1, 7, 'neither'), Interval(3, 5, 'both'), Interval(3, 5, 'both')),        
            # (0, 10] * (1, 2] = (1, 2]
            (Interval(0, 10, 'right'), Interval(1, 2, 'right'), Interval(1, 2, 'right')),        

        ]
        for (interval1, interval2, expected) in cases:
            actual1 = interval1 * interval2
            actual2 = interval2 * interval1
            assert actual1 == expected
            assert actual2 == expected

    def test_subtraction(self):
        cases = [
            # (0, 3) - (2, 6) = (0, 2]
            (Interval(0, 3, 'neither'), Interval(2, 6, 'neither'), Interval(0, 2, 'right')),
            # (1, 10) - (5, 6) = [(1,5],[6,10)] 
            (Interval(1, 10, 'neither'), Interval(5, 6, 'neither'), [Interval(1, 5, 'right'), Interval(6, 10, 'left')]),
            # (1, 10) - [5, 6) = [(1,5),[6,10)] 
            (Interval(1, 10, 'neither'), Interval(5, 6, 'left'), [Interval(1, 5, 'neither'), Interval(6, 10, 'left')]),
            # (1, 2) - (2, 4) = (1,2)
            (Interval(1, 2, 'neither'), Interval(2, 4, 'neither'), (Interval(1, 2, 'neither'))),       
        ]
        for (interval1, interval2, expected) in cases:
            actual1 = interval1 - interval2
            assert actual1 == expected
        
        # (1, 2) - (2, 4) - (3, 5) = (1,2)
        assert Interval(1, 2, 'neither') - Interval(2, 4, 'neither') - Interval(3, 5, 'neither') == Interval(1, 2, 'neither')
        # (1, 10) - (5, 6) - (1,2)= [[2,5], [6,10)] 
        assert Interval(1, 10, 'neither') - Interval(5, 6, 'neither') - Interval(1, 2, 'neither') == [Interval(2, 5, 'both'), Interval(6, 10, 'left')]

class TestIntervalArrayArithmetic:
    def test_addition(self):
        cases = [
            # (1, 2] + (2, 3] + (3, 4] = [(1, 4]]
            (IntervalArray([Interval(1, 2, 'right'), Interval(2, 3, 'right'), Interval(3, 4, 'right')]), IntervalArray([Interval(1,4,'right')])),
            # (1, 2] + (2, 3] + (4, 5] = [(1, 3],(4,5]]
            (IntervalArray([Interval(1, 2, 'right'), Interval(2, 3, 'right'), Interval(4, 5, 'right')]), IntervalArray([Interval(1,3,'right'),Interval(4,5,'right')]))        
        ]
        for (interval_arr, expected) in cases:
            actual = interval_arr.sum()
            tm.assert_interval_array_equal(actual, expected)
    
    def test_multiplication(self):
        cases = [
            # (1, 2] * (2, 3] * (3, 4] = [(0, 0)]
            (IntervalArray([Interval(1, 2, 'right'), Interval(2, 3, 'right'), Interval(3, 4, 'right')]), IntervalArray([Interval(0,0,'neither')])),
            # [1, 2] * [2, 3] * [2, 4] = [[2, 2]]
            (IntervalArray([Interval(1, 2, 'both'), Interval(2, 3, 'both'), Interval(4, 5, 'both')]), IntervalArray([Interval(2,2,'both')]))        
        ]
        for (interval_arr, expected) in cases:
            actual = interval_arr.product()
            tm.assert_interval_array_equal(actual, expected)
