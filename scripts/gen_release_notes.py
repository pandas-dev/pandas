from __future__ import print_function
import sys
import json
from pandas.io.common import urlopen
from datetime import datetime


class Milestone(object):

    def __init__(self, title, number):
        self.title = title
        self.number = number

    def __eq__(self, other):
        if isinstance(other, Milestone):
            return self.number == other.number
        return False


class Issue(object):

    def __init__(self, title, labels, number, milestone, body, state):
        self.title = title
        self.labels = set([x['name'] for x in labels])
        self.number = number
        self.milestone = milestone
        self.body = body
        self.closed = state == 'closed'

    def __eq__(self, other):
        if isinstance(other, Issue):
            return self.number == other.number
        return False


def get_issues():
    all_issues = []
    page_number = 1
    while True:
        iss = _get_page(page_number)
        if len(iss) == 0:
            break
        page_number += 1
        all_issues.extend(iss)
    return all_issues


def _get_page(page_number):
    gh_url = ('https://api.github.com/repos/pandas-dev/pandas/issues?'
              'milestone=*&state=closed&assignee=*&page=%d') % page_number
    with urlopen(gh_url) as resp:
        rs = resp.readlines()[0]
    jsondata = json.loads(rs)
    issues = [Issue(x['title'], x['labels'], x['number'],
                    get_milestone(x['milestone']), x['body'], x['state'])
              for x in jsondata]
    return issues


def get_milestone(data):
    if data is None:
        return None
    return Milestone(data['title'], data['number'])


def collate_label(issues, label):
    lines = []
    for x in issues:
        if label in x.labels:
            lines.append('\t- %s(#%d)' % (x.title, x.number))

    return '\n'.join(lines)


def release_notes(milestone):
    issues = get_issues()

    headers = ['New Features', 'Improvements to existing features',
               'API Changes', 'Bug fixes']
    labels = ['New', 'Enhancement', 'API-Change', 'Bug']

    rs = 'pandas %s' % milestone
    rs += '\n' + ('=' * len(rs))
    rs += '\n\n **Release date:** %s' % datetime.today().strftime('%B %d, %Y')
    for i, h in enumerate(headers):
        rs += '\n\n**%s**\n\n' % h
        l = labels[i]
        rs += collate_label(issues, l)

    return rs

if __name__ == '__main__':

    rs = release_notes(sys.argv[1])
    print(rs)
