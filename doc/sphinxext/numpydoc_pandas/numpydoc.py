# PANDAS HACK: Replace attributes param_list by member_list
def alternative_str(self, indent=0, func_role="obj"):
    ns = {
        'signature': self._str_signature(),
        'index': self._str_index(),
        'summary': self._str_summary(),
        'extended_summary': self._str_extended_summary(),
        'parameters': self._str_param_list('Parameters'),
        'returns': self._str_returns('Returns'),
        'yields': self._str_returns('Yields'),
        'other_parameters': self._str_param_list('Other Parameters'),
        'raises': self._str_param_list('Raises'),
        'warns': self._str_param_list('Warns'),
        'warnings': self._str_warnings(),
        'see_also': self._str_see_also(func_role),
        'notes': self._str_section('Notes'),
        'references': self._str_references(),
        'examples': self._str_examples(),
        'attributes': self._str_member_list('Attributes'),
        'methods': self._str_member_list('Methods'),
    }
    ns = dict((k, '\n'.join(v)) for k, v in ns.items())

    rendered = self.template.render(**ns)
    return '\n'.join(self._str_indent(rendered.split('\n'), indent))


from numpydoc.docscrape_sphinx import SphinxDocString
SphinxDocString.__str__ = alternative_str


def setup(app, *args, **kwargs):
    from numpydoc import setup
    app.add_config_value('numpydoc_attributes_as_param_list', True, True)
    return setup(app, *args, **kwargs)

