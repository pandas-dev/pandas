# DO NOT MERGE
# live working test that tests both scenarios:
# pd.read_csv('https://uname:pwd@fqdn:<port>/fname.csv', verify_ssl=False)
# pd.read_csv('https://fqdn:<port>/fname.csv', username='uname', password='pwd', verify_ssl=False)

import pandas as pd

uname='pandasusr'
pwd='pandaspwd'
url = 'https://{}pandastest.mooo.com:5000/'
verify_ssl=False

def get_df(url, uname, pwd, verify_ssl, pd_read_fn, fname):
  furl = url + fname
  kwargs = {}
  if uname or pwd:
    kwargs['auth']=(uname, pwd)
  if verify_ssl is not None:
    kwargs['verify_ssl']=verify_ssl
  print('\n' +furl)
  df = pd_read_fn(furl, **kwargs)
  if type(df) is list: # html
    df = df[0]
  print(df.to_string(index=False))
  print(df.to_json())

fparams = [(pd.read_csv, 'aaa.csv'),
           (pd.read_json, 'jdoc.json'),
           (pd.read_excel, 'ex_doc.xlsx'),
           (pd.read_html, 'html_file.html')]

for pd_read_fn, fname in fparams:
  u = url.format('{}:{}@'.format(uname, pwd))
  get_df( u, None, None, verify_ssl, pd_read_fn, fname)
  u2 = url.format('')
  get_df( u2, uname, pwd, verify_ssl, pd_read_fn, fname)
