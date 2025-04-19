{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessorcallable:: {{ (module.split('.')[1:] + [objname]) | join('.') }}.__call__
