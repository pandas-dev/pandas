from __future__ import division

import os
import sys
from itertools import cycle

from timeit import Timer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from bokeh.mpl import to_bokeh
from numpy.random import randint

from mpltools import style
style.use('ggplot')

class ExistenceBenchmarks(object):


    def time_py_dict(look_for, look_in):
        df_look_for = pd.DataFrame(look_for, columns=['data'])
        dict_look_in = dict(zip(look_in, look_in))

        def time_this():
            result = df_look_for[[x in dict_look_in for x in df_look_for.data]]
            return result.drop_duplicates().sort('data')
            
        return time_this
        
        
    def time_isin_list(look_for, look_in):
        df_look_for = pd.DataFrame(look_for, columns=['data'])
        list_look_in = list(look_in)
        
        def time_this():
            result = df_look_for[df_look_for.data.isin(list_look_in)]
            return result.drop_duplicates().sort('data')
            
        return time_this
        
        
    def time_isin_dict(look_for, look_in):
        df_look_for = pd.DataFrame(look_for, columns=['data'])
        dict_look_in = dict(zip(look_in, look_in))
        
        def time_this():
            result = df_look_for[df_look_for.data.isin(dict_look_in)]
            return result.drop_duplicates().sort('data')
            
        return time_this
        
        
    def time_isin_series(look_for, look_in):
        series_look_in = pd.Series(look_in)
        df_look_for = pd.DataFrame(look_for, columns=['data'])
        
        def time_this():
            result = df_look_for[df_look_for.data.isin(series_look_in)]
            return result.drop_duplicates().sort('data')
            
        return time_this
        
        
    def time_join(look_for, look_in):
        series_look_in = pd.Series(look_in, index=look_in)
        series_look_in.name = 'series_data'
        df_look_for = pd.DataFrame(look_for, columns=['data'], index=look_for)
        
        def time_this():
            result = df_look_for.join(series_look_in, how='inner')
            return result.drop_duplicates()
            
        return time_this
        
        
    def time_join_no_dups(look_for, look_in):
        series_look_in = pd.Series(look_in, index=look_in)
        series_look_in.name = 'series_data'
        df_look_for = pd.DataFrame(look_for, columns=['data'], index=look_for)
        
        def time_this():
            df_look_for.drop_duplicates(inplace=True)
            series_look_in.drop_duplicates(inplace=True)
            result = df_look_for.join(series_look_in, how='inner')
            return result.sort('data')
            
        return time_this
        
        
    def time_query_in(look_for, look_in):
        series_look_in = pd.Series(look_in)
        series_look_in.name = 'data'
        df_look_for = pd.DataFrame(look_for, columns=['data'])
        
        def time_this():
            # series_look_in is not visible to .query unless defined in local function scope.
            s_look_in = series_look_in
            result = df_look_for.query('data in @s_look_in')
            return result.drop_duplicates().sort('data')
    
        return time_this
        
    
def run_bench(to_time, repeat, look_sets, x_axis, linestyle='-'):
    func_results = []
    markers = cycle(['o', 's', '+', '^', 'v', 'x', 'D', '*'])
    
    for time_func_name in to_time:
        marker=markers.next()
        colors = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
        for set_name, look_set in look_sets:
            color=colors.next()
            plot_results = []
            for look_for, look_in in look_set:
                func = ExistenceBenchmarks.__dict__[time_func_name](look_for, look_in)
                result = func()
                t = Timer(func)
                elapsed = t.timeit(number=repeat) / repeat
                name = time_func_name.replace('time_', '') + ' ' + set_name + ' (%.1f%%)' % ((len(result) / len(look_for)) * 100)
                func_results.append((name, look_for, look_in, elapsed))
                plot_results.append(elapsed)
            plt.plot(x_axis, plot_results, marker=marker, color=color, label=name, linestyle=linestyle)
            
            
def test_timed(to_time):
    look_for = randint(0, 10000, 5000)
    look_in = randint(5000, 15000, 5000)
    
    first_result = ExistenceBenchmarks.__dict__[to_time[0]](look_for, look_in)()

    for time_func_name in to_time[1:]:
        func = ExistenceBenchmarks.__dict__[time_func_name](look_for, look_in)
        result = func()
        if np.array_equal(first_result['data'].values, result['data'].values):
            pass
        else:
            raise AssertionError("%s and %s have unmatched output." % (to_time[0], time_func_name))

        
if __name__ == '__main__':

    pandas_dir = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    static_path = os.path.join(pandas_dir, 'doc', 'source', '_static')
    join_path = lambda p: os.path.join(static_path, p)
    
    to_time = [key for key in ExistenceBenchmarks.__dict__ if key.startswith('time_')]
        
    
    if len(sys.argv) != 2:
        print 'usage: <--test, --run>'
        print '\t--test : Ensure that all timed functions are returning identical output.'
        print '\t--run  : Generate plots for all timed functions.'
        sys.exit()
    
    if sys.argv[1] == '--test':
        test_timed(to_time)
        
    elif sys.argv[1] == '--run':
        test_timed(to_time)
        
        def save_plot(filename, subtitle):
            fname = join_path(filename)
            plt.axes().set_xscale('log')
            x1,x2,y1,y2 = plt.axis()
            # plt.axis((x1, x2, 0, y_limit))
            plt.legend(loc=2, prop={'size':8})
            plt.title('Existence Comparisons%s' % subtitle)
            plt.xlabel('% Overlap of X Elements')
            plt.ylabel('Time(s)')
            plt.savefig(fname)
            plt.clf()
        
        def unordered(exp_range, repeat):   
            rng = [2**x for x in exp_range]
                
            # 25% overlap
            look_set_25 = \
                [(randint(0, 100*i, 50*i), randint(75*i, 175*i, 50*i)) for i in rng]
                
            look_set_50 = \
                [(randint(0, 100*i, 50*i), randint(50*i, 150*i, 50*i)) for i in rng]
                
            look_set_75 = \
                [(randint(0, 100*i, 50*i), randint(25*i, 125*i, 50*i)) for i in rng]
                
            look_set_100 = \
                [(randint(0, 100*i, 50*i), randint(0*i, 100*i, 50*i)) for i in rng]
                
            look_sets = []
            look_sets.append(('25% overlap',  look_set_25))
            look_sets.append(('50% overlap',  look_set_50))
            look_sets.append(('75% overlap',  look_set_75))
            look_sets.append(('100% overlap', look_set_100))
            
            x_axis = [100*i for i in rng]
            run_bench(to_time, 10, look_sets, x_axis, linestyle='-')
        
        
        def from_ordered(exp_range, repeat):
            rng = [2**x for x in exp_range]
                
            # 25% overlap
            look_set_25 = \
                [(sorted(randint(0, 100*i, 50*i)), randint(75*i, 175*i, 50*i)) for i in rng]
                
            look_set_50 = \
                [(sorted(randint(0, 100*i, 50*i)), randint(50*i, 150*i, 50*i)) for i in rng]
                
            look_set_75 = \
                [(sorted(randint(0, 100*i, 50*i)), randint(25*i, 125*i, 50*i)) for i in rng]
                
            look_set_100 = \
                [(sorted(randint(0, 100*i, 50*i)), randint(0*i, 100*i, 50*i)) for i in rng]
                
            look_sets = []
            look_sets.append(('25% overlap, for-ordered',  look_set_25))
            look_sets.append(('50% overlap, for-ordered',  look_set_50))
            look_sets.append(('75% overlap, for-ordered',  look_set_75))
            look_sets.append(('100% overlap, for-ordered', look_set_100))
            
            x_axis = [100*i for i in rng]
            run_bench(to_time, 10, look_sets, x_axis, linestyle='-.')
        

        def both_ordered(exp_range, repeat):
            rng = [2**x for x in exp_range]
                
            # 25% overlap
            look_set_25 = \
                [(sorted(randint(0, 100*i, 50*i)), sorted(randint(75*i, 175*i, 50*i))) for i in rng]
                
            look_set_50 = \
                [(sorted(randint(0, 100*i, 50*i)), sorted(randint(50*i, 150*i, 50*i))) for i in rng]
                
            look_set_75 = \
                [(sorted(randint(0, 100*i, 50*i)), sorted(randint(25*i, 125*i, 50*i))) for i in rng]
                
            look_set_100 = \
                [(sorted(randint(0, 100*i, 50*i)), sorted(randint(0*i, 100*i, 50*i))) for i in rng]
                
            look_sets = []
            look_sets.append(('25% overlap, both-ordered',  look_set_25))
            look_sets.append(('50% overlap, both-ordered',  look_set_50))
            look_sets.append(('75% overlap, both-ordered',  look_set_75))
            look_sets.append(('100% overlap, both-ordered', look_set_100))
            
            x_axis = [100*i for i in rng]
            run_bench(to_time, repeat, look_sets, x_axis, linestyle=':')
            
        
        plt.figure(figsize=(32, 24))
        unordered(range(1, 10), 10)
        from_ordered(range(1, 10), 10)
        both_ordered(range(1, 10), 10)
        save_plot('existence-perf-small.png', ': Small')
        
        plt.figure(figsize=(32, 24))
        unordered(range(10, 15), 3)
        from_ordered(range(10, 15), 3)
        both_ordered(range(10, 15), 3)
        save_plot('existence-perf-large.png', ': Large')
        
        plt.figure(figsize=(16, 12))
        unordered(range(1, 10), 10)
        save_plot('existence-perf-unordered-small.png', ': Unordered Small')
        
        plt.figure(figsize=(16, 12))
        from_ordered(range(1, 10), 10)
        save_plot('existence-perf-from-ordered-small.png', ': From-Ordered Small')
        
        plt.figure(figsize=(16, 12))
        both_ordered(range(1, 10), 10)
        save_plot('existence-perf-both-ordered-small.png', ': Both-Ordered Small')
        
        plt.figure(figsize=(16, 12))
        unordered(range(10, 15), 3)
        save_plot('existence-perf-unordered-large.png', ': Unordered Large')
        
        plt.figure(figsize=(16, 12))
        from_ordered(range(10, 15), 3)
        save_plot('existence-perf-from-ordered-large.png', ': From-Ordered Large')
        
        plt.figure(figsize=(16, 12))
        both_ordered(range(10, 15), 3)
        save_plot('existence-perf-both-ordered-large.png', ': Both-Ordered Large')