from numpy import NaN
from pandas import compat
import numpy as np
import tables
import random

from pandas.core.api import Series, DataFrame, isnull, notnull
from pandas.core.series import remove_na
from pandas.compat import zip

def _data_generate():
	frame1 = random.sample(range(1,100), 1000)
	frame2 = random.sample(range(1,100), 1000)

def _data_MMSTD():
	_data_generate()
	med1 = median(frame1)
	med2 = median(frame2)
	mean1 = mean(frame1)
	mean2 = mean(frame2)
	std1 = stdev(frame1)
	std2 = stdev(frame2)

	text1 = 'Data set 1 has a median of ' + med1 + ', mean of ' + mean1 + ' and standard deviation of ' + std1
	text2 = 'Data set 2 has a median of ' + med2 + ', mean of ' + mean2 + ' and standard deviation of ' + std2

	return text1 + ', and ' + text2


def _data_spread_description():
	_data_generate()	
	_data_MMSTD()

	text3 = ''
	text4 = ''
	text5 = ''

	if(std1 > std2):
	text3 = 'Data set 1 is more spread out than data set 2.'
	else:
	text3 = 'Data set 2 is more spread out than data set 2.'
	
	if(med1 > mean1):
	text4 = 'Data set 1 is left tailed.'
	elif(med1 < mean1):
	text4 = 'Data set 1 is right tailed.'
	else():
	text4 = 'Data set 1 is approximately normally distributed."
	
	if(med2 > mean2):
	text5 = 'Data set 2 is left tailed.'
	elif(med2 < mean2):
	text5 = 'Data set 2 is right tailed.'
	else():
	text5 = 'Data set 2 is approximately normally distributed."

	return text3 + text4 + text5


def _data_stats_tabular():
	_data_generate()
	_data_MMSTD()
	_data_spread_description()
	
	t1 = newTable("Table1")
	t1 = newTable("MyTable",2,4)
	t2 = newTable("Table2")
	t2 = newTable("MyTable1",2,4)

	t1.setText(1, 1, 'Mean')
	t1.setText(2, 1, 'Median')
	t1.setText(3, 1, 'Standard Deviation')
	t1.setText(4, 1, 'Spread')
	t1.setText(2, 1, mean1)
	t1.setText(2, 2, med1)
	t1.setText(2, 3, std1)

	if(text4.find("left tailed"):
	t1.setText(2, 4, 'Left tailed')
	elif(text4.find("right tailed"):
	t1.setText(2, 4, 'Right tailed')
	else():
	t1.setText(2, 4, 'Approximately normal')

	t2.setText(1, 1, 'Mean')
	t2.setText(2, 1, 'Median')
	t2.setText(3, 1, 'Standard Deviation')
	t2.setText(4, 1, 'Spread')
	t2.setText(2, 1, mean2)
	t2.setText(2, 2, med2)
	t2.setText(2, 3, std2)

	if(text5.find("left tailed"):
	t2.setText(2, 4, 'Left tailed')
	elif(text5.find("right tailed"):
	t2.setText(2, 4, 'Right tailed')
	else():
	t2.setText(2, 4, 'Approximately normal')

	print t1
	print t2
	