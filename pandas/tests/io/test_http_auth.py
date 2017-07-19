import pytest
import pandas as pd


def gen_http_auth_ssl_test_cases():
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
        for host, verify_ssl in [('pandas-unittest.site11.com', False, ),
                                 ('handsome-equator.000webhostapp.com', True)
                                 ]:
            for sub_path, is_auth in [('no_auth/', False),
                                      ('basic_auth/', True)]:
                pre_ports = [('http', ''),
                             ('https', '')]
                for pre, port in pre_ports:
                    test_cases.append(
                        [host, verify_ssl, pre, port, sub_path, is_auth])
        return test_cases

    def gen_base_url(pre, auth_prefix, host, port, su_pa):
        return '{}://{}{}{}/{}'.format(pre, auth_prefix, host, port, su_pa)
    tc2 = []
    uname = 'pandasusr'
    pwd = 'pandaspwd'
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


@pytest.mark.parametrize('url, uname, pwd, verify_ssl',
                         gen_http_auth_ssl_test_cases())
def test_http_auth_ssl(url, uname, pwd, verify_ssl):

    def get_df(url, uname, pwd, verify_ssl, pd_read_fn, fname):
        furl = url + fname
        kwargs = {}
        if uname or pwd:
            kwargs['auth'] = (uname, pwd)
        if verify_ssl is not None:
            kwargs['verify_ssl'] = verify_ssl
        msg = '{0: <90} -- auth:[{1: <10}/{2: <10}] v:[{3: <5}]'.format(
              url, str(uname), str(pwd), str(verify_ssl))
        tcsv = 'animal,bird\ndog,pigeon\ncat,emu\n'
        j = '{"animal":{"0":"dog","1":"cat"},"bird":{"0":"pigeon","1":"emu"}}'
        try:
            df = pd_read_fn(furl, **kwargs)
            if type(df) is list:  # html
                df = df[0]
            smatch = str(df.to_csv(index=False)) == tcsv
            jmatch = str(df.to_json()) == j
            res = 'Json : {}  -- String: {}'.format(jmatch, smatch)
            if not jmatch or not smatch:
                raise Exception(' ** ERROR:' + res)
            else:
                res += ' OK'
            print(msg + ' ' + res)
        except Exception as ex:
            print(msg + ' ' + str(ex))
            raise ex
        return True

    for pd_read_fn, fname in [(pd.read_csv, 'aaa.csv'),
                              (pd.read_json, 'jdoc.json'),
                              (pd.read_excel, 'ex_doc.xlsx'),
                              (pd.read_html, 'html_file.html')
                              ]:
        assert get_df(url, uname, pwd, verify_ssl, pd_read_fn, fname)
