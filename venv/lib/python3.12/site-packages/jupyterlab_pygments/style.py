# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from pygments.style import Style
from pygments.token import (
    Comment, Error, Generic, Keyword, Literal, Name, Number, Operator, Other,
    Punctuation, String, Text, Whitespace)


class JupyterStyle(Style):
    """
    A pygments style using JupyterLab CSS variables.

    The goal is to mimick JupyterLab's codemirror theme.

    Known impossibilities:

    - With pygments, the dot in `foo.bar` is considered an Operator (class: 'o'),
      while in codemirror, it is bare text.
    - With pygments, in both `from foo import bar`, and `foo.bar`, "bar" is
      considered a Name (class: 'n'), while in coremirror, the latter is a property.

Available CSS variables are

  --jp-mirror-editor-keyword-color
  --jp-mirror-editor-atom-color
  --jp-mirror-editor-number-color
  --jp-mirror-editor-def-color
  --jp-mirror-editor-variable-color
  --jp-mirror-editor-variable-2-color
  --jp-mirror-editor-variable-3-color
  --jp-mirror-editor-punctuation-color
  --jp-mirror-editor-property-color
  --jp-mirror-editor-operator-color
  --jp-mirror-editor-comment-color
  --jp-mirror-editor-string-color
  --jp-mirror-editor-string-2-color
  --jp-mirror-editor-meta-color
  --jp-mirror-editor-qualifier-color
  --jp-mirror-editor-builtin-color
  --jp-mirror-editor-bracket-color
  --jp-mirror-editor-tag-color
  --jp-mirror-editor-attribute-color
  --jp-mirror-editor-header-color
  --jp-mirror-editor-quote-color
  --jp-mirror-editor-link-color
  --jp-mirror-editor-error-color
    """

    default_style = ''
    background_color = 'var(--jp-cell-editor-background)'
    highlight_color = 'var(--jp-cell-editor-active-background)'

    styles = {
        Text:                      'var(--jp-mirror-editor-variable-color)',        # no class
        Whitespace:                '',                                              # class: 'w'
        Error:                     'var(--jp-mirror-editor-error-color)',           # class: 'err'
        Other:                     '',                                              # class: 'x'

        Comment:                   'italic var(--jp-mirror-editor-comment-color)',  # class: 'c'
        #Comment.Multiline:         '',                                             # class: 'cm'
        #Comment.Preproc:           '',                                             # class: 'cp'
        #Comment.Single:            '',                                             # class: 'c1'
        #Comment.Special:           '',                                             # class: 'cs'

        Keyword:                   'bold var(--jp-mirror-editor-keyword-color)',    # class: 'k'
        #Keyword.Constant:          '',                                             # class: 'kc'
        #Keyword.Declaration:       '',                                             # class: 'kd'
        #Keyword.Namespace:         '',                                             # class: 'kn'
        #Keyword.Pseudo:            '',                                             # class: 'kp'
        #Keyword.Reserved:          '',                                             # class: 'kr'
        #Keyword.Type:              '',                                             # class: 'kt'

        Operator:                  'bold var(--jp-mirror-editor-operator-color)',   # class: 'o'
        Operator.Word:             '',                                              # class: 'ow'

        Literal:                   '',                                              # class: 'l'
        Literal.Date:              '',                                              # class: 'ld'

        String:                    'var(--jp-mirror-editor-string-color)',
        #String.Backtick:           '',                                             # class: 'sb'
        #String.Char:               '',                                             # class: 'sc'
        #String.Doc:                '',                                             # class: 'sd'
        #String.Double:             '',                                             # class: 's2'
        #String.Escape:             '',                                             # class: 'se'
        #String.Heredoc:            '',                                             # class: 'sh'
        #String.Interpol:           '',                                             # class: 'si'
        #String.Other:              '',                                             # class: 'sx'
        #String.Regex:              '',                                             # class: 'sr'
        #String.Single:             '',                                             # class: 's1'
        #String.Symbol:             '',                                             # class: 'ss'

        Number:                    'var(--jp-mirror-editor-number-color)',          # class: 'm'
        #Number.Float:              '',                                             # class: 'mf'
        #Number.Hex:                '',                                             # class: 'mh'
        #Number.Integer:            '',                                             # class: 'mi'
        #Number.Integer.Long:       '',                                             # class: 'il'
        #Number.Oct:                '',                                             # class: 'mo'

        Name:                      '',                                              # class: 'n'
        #Name.Attribute:            '',                                             # class: 'na'
        #Name.Builtin:              '',                                             # class: 'nb'
        #Name.Builtin.Pseudo:       '',                                             # class: 'bp'
        #Name.Class:                '',                                             # class: 'nc'
        #Name.Constant:             '',                                             # class: 'no'
        #Name.Decorator:            '',                                             # class: 'nd'
        #Name.Entity:               '',                                             # class: 'ni'
        #Name.Exception:            '',                                             # class: 'ne'
        #Name.Function:             '',                                             # class: 'nf'
        #Name.Property:             '',                                             # class  'py'
        #Name.Label:                '',                                             # class: 'nl'
        #Name.Namespace:            '',                                             # class: 'nn'
        #Name.Other:                '',                                             # class: 'nx'
        #Name.Tag:                  '',                                             # class: 'nt'
        #Name.Variable:             '',                                             # class: 'nv'
        #Name.Variable.Class:       '',                                             # class: 'vc'
        #Name.Variable.Global:      '',                                             # class: 'vg'
        #Name.Variable.Instance:    '',                                             # class: 'vi'

        Generic:                   '',                                              # class: 'g'
        #Generic.Deleted:           '',                                             # class: 'gd',
        #Generic.Emph:              'italic',                                       # class: 'ge'
        #Generic.Error:             '',                                             # class: 'gr'
        #Generic.Heading:           '',                                             # class: 'gh'
        #Generic.Inserted:          '',                                             # class: 'gi'
        #Generic.Output:            '',                                             # class: 'go'
        #Generic.Prompt:            '',                                             # class: 'gp'
        #Generic.Strong:            '',                                             # class: 'gs'
        #Generic.Subheading:        '',                                             # class: 'gu'
        #Generic.Traceback:         '',                                             # class: 'gt'

        Punctuation:               'var(--jp-mirror-editor-punctuation-color)'       # class: 'p'
    }
