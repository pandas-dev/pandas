from docutils import nodes
from docutils.parsers.rst import Directive

from announce import build_components


class ContributorsDirective(Directive):
    required_arguments = 1
    name = 'contributors'

    def run(self):
        components = build_components(self.arguments[0])

        message = nodes.paragraph()
        message += nodes.Text(components['author_message'])

        listnode = nodes.bullet_list()

        for author in components['authors']:
            para = nodes.paragraph()
            para += nodes.Text(author)
            listnode += nodes.list_item('', para)

        return [message, listnode]


def setup(app):
    app.add_directive('contributors', ContributorsDirective)

    return {'version': '0.1'}
