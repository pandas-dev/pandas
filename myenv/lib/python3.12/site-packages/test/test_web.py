# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import re
import shutil
import time
import urllib.parse
from os.path import join, abspath, dirname

import pytest

from asv import config, util

try:
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.webdriver import ActionChains
    from selenium.common.exceptions import NoSuchElementException, StaleElementReferenceException
    from selenium.webdriver.common.by import By
except ImportError:
    pass

from . import tools
from .tools import get_with_retry, WAIT_TIME, WIN


def _rebuild_basic_html(basedir):
    local = abspath(dirname(__file__))
    cwd = os.getcwd()

    if os.path.isdir(basedir):
        html_dir = join(basedir, 'html')
        dvcs = tools.Git(join(basedir, 'repo'))
        return html_dir, dvcs

    os.makedirs(basedir)
    os.chdir(basedir)
    try:
        machine_file = join(basedir, 'asv-machine.json')

        shutil.copyfile(join(local, 'asv-machine.json'),
                        machine_file)

        values = [[x] * 2 for x in [0, 0, 0, 0, 0,
                                    1, 1, 1, 1, 1,
                                    3, 3, 3, 3, 3,
                                    2, 2, 2, 2, 2]]
        dvcs = tools.generate_test_repo(basedir, values)
        first_tested_commit_hash = dvcs.get_hash(f'{util.git_default_branch()}~14')

        repo_path = dvcs.path
        shutil.move(repo_path, join(basedir, 'repo'))
        dvcs = tools.Git(join(basedir, 'repo'))

        conf = config.Config.from_json({
            'env_dir': join(basedir, 'env'),
            'benchmark_dir': join(local, 'benchmark'),
            'results_dir': join(basedir, 'results_workflow'),
            'html_dir': join(basedir, 'html'),
            'repo': join(basedir, 'repo'),
            'dvcs': 'git',
            'project': 'asv',
            'matrix': {"env": {"SOME_TEST_VAR": ["1"]}},
            'regressions_first_commits': {
                '.*': first_tested_commit_hash
            },
        })

        if WIN:
            # Tell conda to not use hardlinks: on Windows it's not possible
            # to delete hard links to files in use, which causes problem when
            # trying to cleanup environments during this test (since the
            # same cache directory may get reused).
            conf.matrix["env"]["CONDA_ALWAYS_COPY"] = ["True"]

        tools.run_asv_with_conf(conf, 'run', 'ALL',
                                '--show-stderr', '--quick',
                                '--bench=params_examples[a-z0-9_.]*track_',
                                _machine_file=machine_file)

        # Swap CPU info and obtain some results
        info = util.load_json(machine_file, api_version=1)

        # Put in parameter values that need quoting in file names
        info['orangutan']['cpu'] = 'Not /really/ <fast>'
        info['orangutan']['ram'] = '?'
        info['orangutan']['NUL'] = ''

        util.write_json(machine_file, info, api_version=1)

        tools.run_asv_with_conf(conf, 'run', f'{util.git_default_branch()}~10..', '--steps=3',
                                '--show-stderr', '--quick',
                                '--bench=params_examples[a-z0-9_.]*track_',
                                _machine_file=machine_file)

        # Output
        tools.run_asv_with_conf(conf, 'publish')

        shutil.rmtree(join(basedir, 'env'))
    finally:
        os.chdir(cwd)

    return conf.html_dir, dvcs


@pytest.mark.flaky(reruns=1, reruns_delay=5)
def test_web_summarygrid(browser, basic_html):
    html_dir, dvcs = basic_html

    ignore_exc = (NoSuchElementException, StaleElementReferenceException)

    with tools.preview(html_dir) as base_url:
        get_with_retry(browser, base_url)

        WebDriverWait(browser, WAIT_TIME).until(EC.title_is(
            'airspeed velocity of an unladen asv'))

        # Verify benchmark names are displayed as expected
        for href, expected in (
            ('#subdir.time_subdir.time_foo', u'time_subdir.time_foo'),
            ('#params_examples.ParamSuite.track_value', u'ParamSuite.track_value'),
            ('#custom.time_function', u'My Custom Function'),
            ('#named.track_custom_pretty_name', u'this.is/the.answer'),
        ):
            item = browser.find_element(By.XPATH,
                                        f"//a[@href='{href}']/div[@class='benchmark-text']")
            assert item.text == expected

        # Open a graph display, scroll to item and click
        item = browser.find_element(By.LINK_TEXT, 'track_param')

        y = item.location['y']
        browser.execute_script(f'window.scrollTo(0, {y - 200})')

        item.click()

        # Verify there's a plot of some sort
        browser.find_element(By.CSS_SELECTOR, 'canvas.flot-base')

        # Click a parameterized test button, which should toggle the button
        param_button = browser.find_element(By.LINK_TEXT, 'benchmark.params_examples.ClassOne')
        assert 'active' in param_button.get_attribute('class').split()
        param_button.click()

        def check(*args):
            param_button = browser.find_element(By.LINK_TEXT, 'benchmark.params_examples.ClassOne')
            return 'active' not in param_button.get_attribute('class').split()
        WebDriverWait(browser, WAIT_TIME, ignored_exceptions=ignore_exc).until(check)

        # Check there's no error popup; needs an explicit wait because
        # there is no event that occurs on successful load that
        # doesn't also occur on a failed load
        time.sleep(1.0)
        error_box = browser.find_element(By.ID, 'error-message')
        assert not error_box.is_displayed()


@pytest.mark.flaky(reruns=1, reruns_delay=5)
def test_web_regressions(browser, basic_html):
    html_dir, dvcs = basic_html

    bad_commit_hash = dvcs.get_hash(f'{util.git_default_branch()}~9')

    ignore_exc = (NoSuchElementException, StaleElementReferenceException)

    browser.set_window_size(1200, 900)

    with tools.preview(html_dir) as base_url:
        get_with_retry(browser, base_url)

        regressions_btn = browser.find_element(By.LINK_TEXT, 'Regressions')
        regressions_btn.click()

        # Wait for element to appear in the table
        WebDriverWait(browser, WAIT_TIME).until(EC.text_to_be_present_in_element(
            ('xpath', '//table[1]/tbody/tr[2]/td[1]'), 'params_examples.track_find_test'
        ))

        # Check that the expected links appear in the table
        regression_1 = browser.find_element(By.LINK_TEXT, 'params_examples.track_find_test(1)')
        browser.find_element(By.LINK_TEXT, 'params_examples.track_find_test(2)')
        browser.find_element(By.LINK_TEXT, bad_commit_hash[:8])

        href = regression_1.get_attribute('href')
        assert '/#params_examples.track_find_test?' in href
        assert 'commits=' in href

        # Sort the tables vs. benchmark name (PhantomJS doesn't allow doing it via actionchains)
        browser.execute_script("$('thead th').eq(0).stupidsort('asc')")
        WebDriverWait(browser, WAIT_TIME).until(EC.text_to_be_present_in_element(
            ('xpath', '//table[1]/tbody/tr[1]/td[1]'), 'params_examples.track_find_test(1)'
        ))

        # Check the contents of the table
        table_rows = browser.find_elements(By.XPATH, '//table[1]/tbody/tr')
        assert len(table_rows) == 2
        cols1 = [td.text for td in table_rows[0].find_elements(By.XPATH, 'td')]
        cols2 = [td.text for td in table_rows[1].find_elements(By.XPATH, 'td')]

        assert cols1[0] == 'params_examples.track_find_test(1)'
        assert cols2[0] == 'params_examples.track_find_test(2)'

        assert re.match(r'^\d\d\d\d-\d\d-\d\d \d\d:\d\d$', cols1[1])
        assert re.match(r'^\d\d\d\d-\d\d-\d\d \d\d:\d\d$', cols2[1])

        assert cols1[2:] == [bad_commit_hash[:8], '2.00x', '1.00', '2.00', 'Ignore']
        assert cols2[2:] == [bad_commit_hash[:8], '2.00x', '1.00', '2.00', 'Ignore']

        # Check that the ignore buttons work as expected
        buttons = [button for button in browser.find_elements(By.XPATH, '//button')
                   if button.text == 'Ignore']
        buttons[0].click()

        # The button should disappear, together with the link
        WebDriverWait(browser, WAIT_TIME).until_not(EC.visibility_of(buttons[0]))
        WebDriverWait(browser, WAIT_TIME).until_not(EC.visibility_of(regression_1))

        table_rows = browser.find_elements(By.XPATH, '//table[1]/tbody/tr')
        assert len(table_rows) == 1

        # There's a second button for showing the links, clicking
        # which makes the elements reappear
        show_button = [button for button in browser.find_elements(By.XPATH, '//button')
                       if button.text == 'Show ignored regressions...'][0]
        show_button.click()

        regression_1 = browser.find_element(By.LINK_TEXT, 'params_examples.track_find_test(1)')
        WebDriverWait(browser, WAIT_TIME).until(EC.visibility_of(regression_1))

        table_rows = browser.find_elements(By.XPATH, '//table[2]/tbody/tr')
        assert len(table_rows) == 1

        # There's a config sample element
        pre_div = browser.find_element(By.XPATH, '//pre')
        assert "params_examples\\\\.track_find_test\\\\(1\\\\)" in pre_div.text

        # There's an unignore button that moves the element back to the main table
        unignore_button = [button for button in browser.find_elements(By.XPATH, '//button')
                           if button.text == 'Unignore'][0]
        unignore_button.click()

        # wait until the table has two rows
        browser.find_elements(By.XPATH, '//table[1]/tbody/tr[2]')

        table_rows = browser.find_elements(By.XPATH, '//table[1]/tbody/tr')
        assert len(table_rows) == 2

        # Check that a plot of some sort appears on mouseover.  The
        # page needs to be scrolled first so that the mouseover popup
        # has enough space to appear.
        regression_1 = browser.find_element(By.LINK_TEXT, 'params_examples.track_find_test(1)')

        y = regression_1.location['y']
        browser.execute_script(f'window.scrollTo(0, {y - 200})')

        chain = ActionChains(browser)
        chain.move_to_element(regression_1)
        chain.perform()

        browser.find_element(By.CSS_SELECTOR, 'div.popover-content')
        browser.find_element(By.CSS_SELECTOR, 'canvas.flot-base')

        # Check group/ungroup button functionality
        group_button, = [button for button in browser.find_elements(By.XPATH, '//button')
                         if button.text == "Group regressions"]
        group_button.click()

        def check(*args):
            columns = browser.find_element(By.XPATH, '//table/thead/tr[1]').text
            return columns == 'Benchmark Last date Commits Factor Best Current'

        WebDriverWait(browser, WAIT_TIME, ignored_exceptions=ignore_exc).until(check)

        ungroup_button, = [button for button in browser.find_elements(By.XPATH, '//button')
                           if button.text == "Ungroup regressions"]
        ungroup_button.click()

        def check(*args):
            columns = browser.find_element(By.XPATH, '//table/thead/tr[1]').text
            return columns == 'Benchmark Date Commit Factor Before Best after'

        WebDriverWait(browser, WAIT_TIME, ignored_exceptions=ignore_exc).until(check)


@pytest.mark.flaky(reruns=1, reruns_delay=5)
def test_web_summarylist(browser, basic_html):
    ignore_exc = (NoSuchElementException, StaleElementReferenceException)

    html_dir, dvcs = basic_html

    last_change_hash = dvcs.get_hash(f'{util.git_default_branch()}~4')

    browser.set_window_size(1200, 900)

    with tools.preview(html_dir) as base_url:
        get_with_retry(browser, base_url)

        summarylist_btn = browser.find_element(By.LINK_TEXT, 'Benchmark list')
        summarylist_btn.click()

        # Check text content in the table
        base_link = browser.find_element(By.LINK_TEXT, 'params_examples.track_find_test')
        cur_row = base_link.find_element(By.XPATH, '../..')
        m = re.match('params_examples.track_find_test \\([12]\\) 2.00 \u221233.3% \\(-1.00\\).*' +
                     last_change_hash[:8],
                     cur_row.text)
        assert m, cur_row.text

        # Check units in row
        base_link2 = browser.find_element(By.LINK_TEXT, 'params_examples.track_bytes')
        cur_row2 = base_link2.find_element(By.XPATH, '../..')
        m = re.match(r'params_examples.track_bytes\s*1.000M', cur_row2.text)
        assert m, cur_row2.text

        # Check link
        base_href, qs = urllib.parse.splitquery(base_link.get_attribute('href'))
        base_url, tag = urllib.parse.splittag(base_href)
        assert urllib.parse.parse_qs(qs) == {'ram': ['128GB'], 'cpu': ['Blazingly fast'],
                                             'NUL': ['[none]']}
        assert tag == 'params_examples.track_find_test'

        # Change table sort (sorting is async, so needs waits)
        sort_th = browser.find_element(By.XPATH, '//th[text()="Recent change"]')
        sort_th.click()
        WebDriverWait(browser, WAIT_TIME).until(
            EC.text_to_be_present_in_element(('xpath', '//tbody/tr[1]'),
                                             'params_examples.track_find_test'))

        # Try to click cpu selector link in the panel
        cpu_select = browser.find_element(By.LINK_TEXT, 'Not /really/ <fast>')
        cpu_select.click()

        # For the other CPU, there is no recent change recorded, only
        # the latest result is available
        def check(*args):
            links = browser.find_elements(By.LINK_TEXT, 'params_examples.track_find_test')
            visible_links = [item for item in links if item.is_displayed()]

            row_texts = [link.find_element(By.XPATH, '../..').text
                         for link in visible_links]
            row_texts.sort()

            if len(row_texts) != 2:
                return False

            ok = (re.match(r'^params_examples\.track_find_test \(1\) 2\.00 .*\(-1\.00\).*$',
                           row_texts[0]) and
                  re.match(r'^params_examples\.track_find_test \(2\) 2\.00 .*\(-1\.00\).*$',
                  row_texts[1]))
            return ok

        WebDriverWait(browser, WAIT_TIME, ignored_exceptions=ignore_exc).until(check)
