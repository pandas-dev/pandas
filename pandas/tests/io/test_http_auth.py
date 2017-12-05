import pytest
import pandas as pd
import pandas.util.testing as tm

if pd.io.common._REQUESTS_INSTALLED:
    from requests.packages.urllib3.exceptions import InsecureRequestWarning
    import requests
    requests_pkg = True
else:
    requests_pkg = False


def _skip_if_no_requests():
    if not requests_pkg:
        pytest.skip('python-requests not installed, skipping')


uname = 'pandasusr'
pwd = 'pandaspwd'
no_auth_path = 'no_auth/'
basic_auth_path = 'basic_auth/'
valid_ssl_url = 'handsome-equator.000webhostapp.com'
invalid_ssl_url = 'pandas-unittest.site11.com'


def gen_http_auth_ssl_test_cases(uname,
                                 pwd,
                                 is_auth,
                                 sub_path):
    """
    Generate list of test case to test for : http/https, username/pwd in url
    or as parameters, self signed ssl certs or trusted ssl certs, no auth
    or basic auth
    """
    def gen_level1_tc():
        test_cases = []
        # The following host doesn't seem to handle urllib but handles
        # python requests package. This is because:
        # 'urlopen' sets header 'Host' : '<fqdn>:<port>' - acceptable RFC7230
        # 'requests' sets header 'Host' : '<fqdn>'
        # so pandas fails on following hosting server (uses some 'apex' server)
        # but pandas succeeds on nginx even if port is non-default.
        for host, verify_ssl in [(invalid_ssl_url, False),
                                 (valid_ssl_url, True)
                                 ]:
            pre_ports = [('http', ''),
                         ('https', '')]
            for pre, port in pre_ports:
                test_cases.append(
                    [host, verify_ssl, pre, port, sub_path, is_auth])
        return test_cases

    def gen_base_url(pre, auth_prefix, host, port, su_pa):
        return '{}://{}{}{}/{}'.format(pre, auth_prefix, host, port, su_pa)
    tc2 = []
    for host, verify_ssl, pre, port, sp, is_auth in gen_level1_tc():
        u = uname if is_auth else None
        p = pwd if is_auth else None
        u_no_uname = gen_base_url(pre, '', host, port, sp)
        u_with_uname = None
        if is_auth:
            auth_prefix = '{}:{}@'.format(u, p) if is_auth else ''
            u_with_uname = gen_base_url(pre, auth_prefix, host, port, sp)
            tc2.append([u_no_uname, u, p, verify_ssl])
            if u_with_uname and u_with_uname != u_no_uname:
                tc2.append([u_with_uname, None, None, verify_ssl])
        else:
            tc2.append([u_no_uname, None, None, verify_ssl])
    return tc2


valid_no_auth = gen_http_auth_ssl_test_cases(uname='',
                                             pwd='',
                                             is_auth=False,
                                             sub_path=no_auth_path)

valid_auth = gen_http_auth_ssl_test_cases(uname=uname,
                                          pwd=pwd,
                                          is_auth=True,
                                          sub_path=basic_auth_path)


@pytest.mark.slow
@pytest.mark.parametrize('url, uname, pwd, verify_ssl',
                         valid_auth)
def test_http_valid_auth(url, uname, pwd, verify_ssl):
    _skip_if_no_requests()
    check_http_auth(url, uname, pwd, verify_ssl)


@pytest.mark.slow
@pytest.mark.parametrize('url, uname, pwd, verify_ssl',
                         valid_no_auth)
def test_http_valid_no_auth(url, uname, pwd, verify_ssl):
    if verify_ssl is False:
        _skip_if_no_requests()
    if (verify_ssl in [True, False] or uname or pwd) and not requests_pkg:
        with pytest.raises(ValueError):
            check_http_auth(url, uname, pwd, verify_ssl)
    else:
        check_http_auth(url, uname, pwd, verify_ssl)


wrong_auth = gen_http_auth_ssl_test_cases(uname='fakepwd',
                                          pwd='fakepwd',
                                          is_auth=True,
                                          sub_path=basic_auth_path)


@pytest.mark.slow
@pytest.mark.parametrize('url, uname, pwd, verify_ssl',
                         wrong_auth)
def test_http_invalid_auth(url, uname, pwd, verify_ssl):
    _skip_if_no_requests()
    with pytest.raises(requests.exceptions.HTTPError):
        check_http_auth(url, uname, pwd, verify_ssl)


blank_uname = gen_http_auth_ssl_test_cases(uname='',
                                           pwd='fakepwd',
                                           is_auth=True,
                                           sub_path=basic_auth_path)

blank_pwd = gen_http_auth_ssl_test_cases(uname='fakepwd',
                                         pwd='',
                                         is_auth=True,
                                         sub_path=basic_auth_path)


def match_csv(df):
    tcsv = 'animal,bird\ndog,pigeon\ncat,emu\n'
    return str(df.to_csv(index=False)) == tcsv


def match_json(df):
    j = '{"animal":{"0":"dog","1":"cat"},"bird":{"0":"pigeon","1":"emu"}}'
    return str(df.to_json()) == j


@tm.network
def check_http_auth(url, uname, pwd, verify_ssl):

    def get_df(url, uname, pwd, verify_ssl, pd_read_fn, fname):
        furl = url + fname
        http_params = {}
        if uname or pwd:
            http_params['auth'] = (uname, pwd)
        if verify_ssl is not None:
            http_params['verify'] = verify_ssl
        msg = '{0: <90} -- auth:[{1: <10}/{2: <10}] v:[{3: <5}]'.format(
              furl, str(uname), str(pwd), str(verify_ssl))
        if verify_ssl or furl.lower().startswith('http://'):
            df = pd_read_fn(furl, http_params=http_params)
        else:
            with tm.assert_produces_warning(InsecureRequestWarning):
                df = pd_read_fn(furl, http_params=http_params)
        if type(df) is list:  # html
            df = df[0]
        smatch = match_csv(df)
        jmatch = match_json(df)
        res = 'Json : {}  -- String: {}'.format(jmatch, smatch)
        if not jmatch or not smatch:
            raise Exception(' ** ERROR:' + res)
        else:
            res += ' OK'
        print(msg + ' ' + res)
        return True

    for pd_read_fn, fname in [(pd.read_csv, 'aaa.csv'),
                              (pd.read_json, 'jdoc.json'),
                              (pd.read_excel, 'ex_doc.xlsx'),
                              (pd.read_html, 'html_file.html')
                              ]:
        assert get_df(url, uname, pwd, verify_ssl, pd_read_fn, fname)
    return


@pytest.mark.parametrize('up',
                         [
                             None,
                             {},
                             {'auth': ('uname', 'pwd')},
                             {'verify': False},
                             {'auth': ('uname', 'pwd'), 'verify': False},
                         ]
                         )
def test_http_params(up):
    _skip_if_no_requests()
    s = pd.io.common.gen_session(http_params=up)
    assert type(s) is requests.sessions.Session
    if up and up.get('auth', None):
        assert s.auth == up.get('auth', None)
    if up and (up.get('verify', True) is False):
        assert s.verify == up.get('verify', True)


def test_pass_session_obj():
    _skip_if_no_requests()
    s = requests.sessions.Session()
    s.auth = ('uname', 'pwd')
    s.verify = False
    t = pd.io.common.gen_session(http_params=s)
    assert s == t
    assert s.auth == t.auth
    assert s.verify == t.verify


def test_skip_requests():
    with pytest.raises(ValueError):
        a = (uname, pwd)
        resp, content_bytes = pd.io.common.fetch_url('http://cnn.com',
                                                     http_params={'auth': a},
                                                     skip_requests=True)
