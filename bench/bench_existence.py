from timeit import Timer
import pandas as pd
import matplotlib.pyplot as plt
import os


class Benchmarks(object):

    def removed_time_py_list(look_for, look_in):
        l = range(look_in)
        df = pd.DataFrame(range(look_for))

        def time_this():
            df[[x in l for x in df.index.values]]
            
        return time_this

    def time_py_dict(look_for, look_in):
        l = range(look_in)
        l_dict = dict(zip(l, l))
        df = pd.DataFrame(range(look_for))

        def time_this():
            df[[x in l_dict for x in df.index.values]]
            
        return time_this
        
        
    def time_isin_list(look_for, look_in):
        l = range(look_in)
        df = pd.DataFrame(range(look_for))
        
        def time_this():
            df[df.index.isin(l)]
            
        return time_this
        
        
    def time_isin_dict(look_for, look_in):
        l = range(look_in)
        l_dict = dict(zip(l, l))
        df = pd.DataFrame(range(look_for))
        
        def time_this():
            df[df.index.isin(l_dict)]
            
        return time_this
        
        
    def time_isin_series(look_for, look_in):
        l = range(look_in)
        l_series = pd.Series(l)
        df = pd.DataFrame(range(look_for))
        
        def time_this():
            df[df.index.isin(l_series.index)]
            
        return time_this
        
        
    def time_join(look_for, look_in):
        l = range(look_in)
        l_series = pd.Series(l)
        l_series.name = 'data'
        df = pd.DataFrame(range(look_for))
        
        def time_this():
            df.join(l_series, how='inner')
            
        return time_this
    
    # Removed. This functionality might be a bug in query('.. == ..').
    # def time_query_eqeq(look_for, look_in):
        # l = range(look_in)
        # s = pd.Series(l)
        # s.name = 'data'
        # df = pd.DataFrame(range(look_for))
        
        # def time_this():
            # l_series = s
            # df.query('index == @l_series')
    
        # return time_this
        
    def time_query_in(look_for, look_in):
        l = range(look_in)
        s = pd.Series(l)
        s.name = 'data'
        df = pd.DataFrame(range(look_for))
        
        def time_this():
            l_series = s
            df.query('index in @l_series')
    
        return time_this
        
    
def run_bench(to_time, repeat, look_in, num_look_for_rows, y_limit, filename):
    func_results = []
    plt.figure()
    
    for time_func_name in to_time:
        plot_results = []
        for look_for in num_look_for_rows:
            func = Benchmarks.__dict__[time_func_name](look_for, look_in)
            t = Timer(func)
            elapsed = t.timeit(number=repeat) / repeat
            name = time_func_name.replace('time_', '')
            func_results.append((name, look_for, look_in, elapsed))
            plot_results.append(elapsed)
        plt.plot(num_look_for_rows, plot_results, label=name)
        
    plt.axes().set_xscale('log')
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1, x2, 0, y_limit))
        
    plt.legend(loc=2, prop={'size':8})
    plt.title('Look in %s Rows' % look_in)
    plt.xlabel('Look For X Rows')
    plt.ylabel('Time(s)')
    plt.savefig(filename)
    plt.clf()
            

if __name__ == '__main__':

    pandas_dir = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    static_path = os.path.join(pandas_dir, 'doc', 'source', '_static')
    
    join = lambda p: os.path.join(static_path, p)
    
    to_time = [key for key in Benchmarks.__dict__ if key.startswith('time_')]
        
    num_look_for_rows = [10 * 2**i for i in range(1, 21)]
        
    filename = join('existence-perf-small.png')
    run_bench(to_time, 10, 5000, num_look_for_rows[0:len(num_look_for_rows)/2], 0.004, filename)
    
    filename = join('existence-perf-large.png')
    run_bench(to_time, 3, 5000000, num_look_for_rows[len(num_look_for_rows)/2:], 10, filename)
    
