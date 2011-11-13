from pandas.core.reshape import melt
import pandas.util.testing as tm

def test_melt():
    df = tm.makeTimeDataFrame()[:10]
    df['id1'] = (df['A'] > 0).astype(int)
    df['id2'] = (df['B'] > 0).astype(int)

    molten1 = melt(df)
    molten2 = melt(df, id_vars=['id1'])
    molten3 = melt(df, id_vars=['id1', 'id2'])

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

