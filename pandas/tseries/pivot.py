import numpy as np

from pandas.core.frame import DataFrame
import pandas.core.nanops as nanops
from pandas.tseries.util import isleapyear
from pandas.tseries.index import date_range

def pivot_annual_h(series, freq=None, dt_index=False):
    """
    Group a series by years, taking leap years into account.

    The output has as many rows as distinct years in the original series,
    and as many columns as the length of a leap year in the units corresponding
    to the original frequency (366 for daily frequency, 366*24 for hourly...).
    The fist column of the output corresponds to Jan. 1st, 00:00:00,
    while the last column corresponds to Dec, 31st, 23:59:59.
    Entries corresponding to Feb. 29th are masked for non-leap years.

    For example, if the initial series has a daily frequency, the 59th column
    of the output always corresponds to Feb. 28th, the 61st column to Mar. 1st,
    and the 60th column is masked for non-leap years.
    With a hourly initial frequency, the (59*24)th column of the output always
    correspond to Feb. 28th 23:00, the (61*24)th column to Mar. 1st, 00:00, and
    the 24 columns between (59*24) and (61*24) are masked.

    If the original frequency is less than daily, the output is equivalent to
    ``series.convert('A', func=None)``.

    Parameters
    ----------
    series : TimeSeries
    freq : string or None, default None

    Returns
    -------
    annual : DataFrame
    
    
    """
    #TODO: test like original pandas and the position of first and last value in arrays
    #TODO: reduce number of hardcoded values scattered all around.   
    index = series.index
    year = index.year
    years = nanops.unique1d(year)    
    
    if freq is not None:
        freq = freq.upper()
    else:
        freq = series.index.freq

    if freq == 'H':
    
        ##basics
        
        #integer value of sum of all hours in a leap hear
        total_hoy_leap = (year_length(series.index.freqstr))
        
        #list of all hours in a leap year
        hoy_leap_list = range(1, (total_hoy_leap + 1 ))
        
        
        #create a array template
        values = np.empty((total_hoy_leap, len(years)), dtype=series.dtype)
        values.fill(np.nan)
        #create a df to receive the resulting data
        dummy_df = DataFrame(values, index=hoy_leap_list, 
                        columns=years)
                        
        ##prepare the index for inserting the values into the result dataframe
        #get offset for leap hours
        #see:
        #http://stackoverflow.com/questions/2004364/increment-numpy-array-with-repeated-indices
        #1994-02-28 23:00:00 -> index 1415
        index_nonleap = np.array(range(0, 8760))
        index_leapshift = np.array(range(1416,8760 ))
        
        index_incl_leap = index_nonleap.copy()
        #shift index by 24 (hours) for leap
        index_incl_leap[index_leapshift]+=24
        
        # select data for the respective year
        for year in years:
            
            #select the data for the respective year
            series_year = series[ series.index.year == year]
            #create a array with the values for the respecive year
            values = (series_year).values
            
            if isleapyear(year):
                dummy_df[year] = values
            else:
                #dummy array to be filled with non-leap values
                dummy_array = np.empty((total_hoy_leap), dtype=series.dtype)
                dummy_array.fill(np.nan)
                
                #fill dummy array with values leaving the leap day
                dummy_array.put(index_incl_leap, values)
                
                dummy_df[year] = dummy_array
                
        res_df = dummy_df
        
        #assign a pseudo datetime index , CAUTION: the year is definitely wrong!
        if dt_index:
            rng = default_rng(freq='H', leap=True)            
            res_df = DataFrame(res_df.values, index=rng, 
                               columns=res_df.columns)
        
        return res_df
        
#TDOO: use pivot_annual for D & M and minute in the same fashion
    if freq == 'D':
        raise NotImplementedError(freq), "use pandas.tseries.util.pivot_annual"        
        
    if freq == 'M':
        raise NotImplementedError(freq), "use pandas.tseries.util.pivot_annual"
    
    else:
        raise NotImplementedError(freq)
        
    
    return res_df
    
    
### timeseries pivoting helper

def last_col2front(df, col_no=1):
    """shifts the last column of a data frame to the front
    
    increase col_no to shift more cols    
    """
    cols = cols = df.columns.tolist()
    #increase index value to 2+ if more columns are to be shifted
    cols = cols[-col_no:] + cols[:-col_no]
    df = df[cols]
    
    return df
    

def extended_info(df, time_cols=True, aggreg=True, aggreg_func=None,
                  datetime_index=False):
    """add extended information to a timeseries pivot
    """

    df_extended = df.copy()
    #perform the following only on the data columns
    cols = df_extended.columns
    #TODO: add standard aggregation
    #TODO: make function be set by argument
    #TODO: is there no a SM describe function?
    #TODO: Maybe use http://pandas.pydata.org/pandas-docs/dev/basics.html#summarizing-data-describe
    if aggreg:
           
        df_extended['mean'] = df_extended[cols].mean(1)
        df_extended['sum'] = df_extended[cols].sum(1)
        df_extended['min'] = df_extended[cols].min(1)
        df_extended['max'] = df_extended[cols].max(1)
        df_extended['std'] = df_extended[cols].std(1)

    #add some metadata
    #TODO: add function to make index a datetime with the argument above using the rng below    
    #TODO: convert the range to lower frequencies and reuse the function.
    rng = default_rng()
    df_extended['doy'] = rng.dayofyear
#    df_extended = last_col2front(df_extended)
    df_extended['month'] = rng.month
#    df_extended = last_col2front(df_extended)
    df_extended['day'] = rng.day
#    df_extended = last_col2front(df_extended)
    df_extended['hour'] = rng.hour + 1
    df_extended = last_col2front(df_extended, col_no=4)
    
    return df_extended
    
###Timeseries convenience / helper functions

                        
def year_length(freq, leap=True):
    """helper function for year length at different frequencies.
    to be expanded
    """

    daysofyear_leap = 366
    daysofyear_nonleap = 365
        
    if freq == 'H':
        if leap:        
            length = 24 * daysofyear_leap
        else:
            length = 24 * daysofyear_nonleap
            
    return length

def default_rng(freq='H', leap=True):
    """create default ranges
    """
    
    if leap:
        total_hoy_leap = (year_length(freq='H'))    
        rng = date_range('1/1/2012', periods=total_hoy_leap, freq='H')
    
    else:
        total_hoy_nonleap = (year_length(freq='H'))    
        rng = date_range('1/1/2011', periods=total_hoy_nonleap, freq='H')        
    
    return rng
