'''
Easy integration of DataFrame into pyqt framework

@author: Jev Kuznetsov
'''
from PyQt4.QtCore import (QAbstractTableModel,Qt,QVariant,QModelIndex, SIGNAL)
from PyQt4.QtGui import (QApplication,QDialog,QVBoxLayout, QTableView, QWidget)

from pandas import DataFrame, Index



class DataFrameModel(QAbstractTableModel):
    ''' data model for a DataFrame class '''
    def __init__(self):
        super(DataFrameModel,self).__init__()
        self.df = DataFrame()
         
    def setDataFrame(self,dataFrame):
        self.df = dataFrame
    
    def signalUpdate(self):
        ''' tell viewers to update their data (this is full update, not efficient)'''
        self.layoutChanged.emit()
              
    #------------- table display functions -----------------     
    def headerData(self,section,orientation,role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return QVariant()
          
        if orientation == Qt.Horizontal:
            try:
                return self.df.columns.tolist()[section]
            except (IndexError, ):
                return QVariant()
        elif orientation == Qt.Vertical:
            try:
                #return self.df.index.tolist()
                return self.df.index.tolist()[section]
            except (IndexError, ):
                return QVariant()
            
    def data(self, index, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return QVariant()
        
        if not index.isValid():
            return QVariant()
        
        return QVariant(str(self.df.ix[index.row(),index.column()]))
        
      
    def rowCount(self, index=QModelIndex()):
        return self.df.shape[0]
    
    def columnCount(self, index=QModelIndex()):
        return self.df.shape[1] 


class DataFrameWidget(QWidget):
    ''' a simple widget for using DataFrames in a gui '''
    def __init__(self,dataFrame, parent=None):
        super(DataFrameWidget,self).__init__(parent)
        
        self.dataModel = DataFrameModel()
        self.dataModel.setDataFrame(dataFrame)
        
        self.dataTable = QTableView()
        self.dataTable.setModel(self.dataModel)
        self.dataModel.signalUpdate()
        
        layout = QVBoxLayout()
        layout.addWidget(self.dataTable)
        self.setLayout(layout)
        
        
    
    def resizeColumnsToContents(self):
        self.dataTable.resizeColumnsToContents()    
        
#-----------------stand alone test code

def testDf():
    ''' creates test dataframe '''
    data = {'int':[1,2,3],'float':[1.5,2.5,3.5],'string':['a','b','c'],'nan':[np.nan,np.nan,np.nan]}
    return DataFrame(data, index=Index(['AAA','BBB','CCC']))[['int','float','string','nan']]


class Form(QDialog):
    def __init__(self,parent=None):
        super(Form,self).__init__(parent)
         
        df = testDf() # make up some data
        widget = DataFrameWidget(df)
        widget.resizeColumnsToContents()
                     
        layout = QVBoxLayout()
        layout.addWidget(widget)
        self.setLayout(layout)
        
if __name__=='__main__':
    import sys
    import numpy as np
    
    app = QApplication(sys.argv)
    form = Form()
    form.show()
    app.exec_()
    

        
        

        