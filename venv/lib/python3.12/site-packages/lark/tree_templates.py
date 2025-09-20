"""This module defines utilities for matching and translation tree templates.

A tree templates is a tree that contains nodes that are template variables.

"""

from typing import Union, Optional, Mapping, Dict, Tuple, Iterator

from lark import Tree, Transformer
from lark.exceptions import MissingVariableError

Branch = Union[Tree[str], str]
TreeOrCode = Union[Tree[str], str]
MatchResult = Dict[str, Tree]
_TEMPLATE_MARKER = '$'


class TemplateConf:
    """Template Configuration

    Allows customization for different uses of Template

    parse() must return a Tree instance.
    """

    def __init__(self, parse=None):
        self._parse = parse

    def test_var(self, var: Union[Tree[str], str]) -> Optional[str]:
        """Given a tree node, if it is a template variable return its name. Otherwise, return None.

        This method may be overridden for customization

        Parameters:
            var: Tree | str - The tree node to test

        """
        if isinstance(var, str):
            return _get_template_name(var)

        if (
            isinstance(var, Tree)
            and var.data == "var"
            and len(var.children) > 0
            and isinstance(var.children[0], str)
        ):
            return _get_template_name(var.children[0])

        return None

    def _get_tree(self, template: TreeOrCode) -> Tree[str]:
        if isinstance(template, str):
            assert self._parse
            template = self._parse(template)

        if not isinstance(template, Tree):
            raise TypeError("template parser must return a Tree instance")

        return template

    def __call__(self, template: Tree[str]) -> 'Template':
        return Template(template, conf=self)

    def _match_tree_template(self, template: TreeOrCode, tree: Branch) -> Optional[MatchResult]:
        """Returns dict of {var: match} if found a match, else None
        """
        template_var = self.test_var(template)
        if template_var:
            if not isinstance(tree, Tree):
                raise TypeError(f"Template variables can only match Tree instances. Not {tree!r}")
            return {template_var: tree}

        if isinstance(template, str):
            if template == tree:
                return {}
            return None

        assert isinstance(template, Tree) and isinstance(tree, Tree), f"template={template} tree={tree}"

        if template.data == tree.data and len(template.children) == len(tree.children):
            res = {}
            for t1, t2 in zip(template.children, tree.children):
                matches = self._match_tree_template(t1, t2)
                if matches is None:
                    return None

                res.update(matches)

            return res

        return None


class _ReplaceVars(Transformer[str, Tree[str]]):
    def __init__(self, conf: TemplateConf, vars: Mapping[str, Tree[str]]) -> None:
        super().__init__()
        self._conf = conf
        self._vars = vars

    def __default__(self, data, children, meta) -> Tree[str]:
        tree = super().__default__(data, children, meta)

        var = self._conf.test_var(tree)
        if var:
            try:
                return self._vars[var]
            except KeyError:
                raise MissingVariableError(f"No mapping for template variable ({var})")
        return tree


class Template:
    """Represents a tree template, tied to a specific configuration

    A tree template is a tree that contains nodes that are template variables.
    Those variables will match any tree.
    (future versions may support annotations on the variables, to allow more complex templates)
    """

    def __init__(self, tree: Tree[str], conf: TemplateConf = TemplateConf()):
        self.conf = conf
        self.tree = conf._get_tree(tree)

    def match(self, tree: TreeOrCode) -> Optional[MatchResult]:
        """Match a tree template to a tree.

        A tree template without variables will only match ``tree`` if it is equal to the template.

        Parameters:
            tree (Tree): The tree to match to the template

        Returns:
            Optional[Dict[str, Tree]]: If match is found, returns a dictionary mapping
                template variable names to their matching tree nodes.
                If no match was found, returns None.
        """
        tree = self.conf._get_tree(tree)
        return self.conf._match_tree_template(self.tree, tree)

    def search(self, tree: TreeOrCode) -> Iterator[Tuple[Tree[str], MatchResult]]:
        """Search for all occurrences of the tree template inside ``tree``.
        """
        tree = self.conf._get_tree(tree)
        for subtree in tree.iter_subtrees():
            res = self.match(subtree)
            if res:
                yield subtree, res

    def apply_vars(self, vars: Mapping[str, Tree[str]]) -> Tree[str]:
        """Apply vars to the template tree
        """
        return _ReplaceVars(self.conf, vars).transform(self.tree)


def translate(t1: Template, t2: Template, tree: TreeOrCode):
    """Search tree and translate each occurrence of t1 into t2.
    """
    tree = t1.conf._get_tree(tree)      # ensure it's a tree, parse if necessary and possible
    for subtree, vars in t1.search(tree):
        res = t2.apply_vars(vars)
        subtree.set(res.data, res.children)
    return tree


class TemplateTranslator:
    """Utility class for translating a collection of patterns
    """

    def __init__(self, translations: Mapping[Template, Template]):
        assert all(isinstance(k, Template) and isinstance(v, Template) for k, v in translations.items())
        self.translations = translations

    def translate(self, tree: Tree[str]):
        for k, v in self.translations.items():
            tree = translate(k, v, tree)
        return tree


def _get_template_name(value: str) -> Optional[str]:
    return value.lstrip(_TEMPLATE_MARKER) if value.startswith(_TEMPLATE_MARKER) else None
