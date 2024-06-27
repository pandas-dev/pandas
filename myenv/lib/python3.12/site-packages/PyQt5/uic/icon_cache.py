#############################################################################
##
## Copyright (c) 2023 Riverbank Computing Limited <info@riverbankcomputing.com>
## 
## This file is part of PyQt5.
## 
## This file may be used under the terms of the GNU General Public License
## version 3.0 as published by the Free Software Foundation and appearing in
## the file LICENSE included in the packaging of this file.  Please review the
## following information to ensure the GNU General Public License version 3.0
## requirements will be met: http://www.gnu.org/copyleft/gpl.html.
## 
## If you do not wish to use this file under the terms of the GPL version 3.0
## then you may purchase a commercial license.  For more information contact
## info@riverbankcomputing.com.
## 
## This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
## WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
##
#############################################################################


import os.path


class IconCache(object):
    """Maintain a cache of icons.  If an icon is used more than once by a GUI
    then ensure that only one copy is created.
    """

    def __init__(self, object_factory, qtgui_module):
        """Initialise the cache."""

        self._object_factory = object_factory
        self._qtgui_module = qtgui_module
        self._base_dir = ''
        self._cache = []

    def set_base_dir(self, base_dir):
        """ Set the base directory to be used for all relative filenames. """

        self._base_dir = base_dir

    def get_icon(self, iconset):
        """Return an icon described by the given iconset tag."""

        # Handle a themed icon.
        theme = iconset.attrib.get('theme')
        if theme is not None:
            return self._object_factory.createQObject("QIcon.fromTheme",
                    'icon', (self._object_factory.asString(theme), ),
                    is_attribute=False)

        # Handle an empty iconset property.
        if iconset.text is None:
            return None

        iset = _IconSet(iconset, self._base_dir)

        try:
            idx = self._cache.index(iset)
        except ValueError:
            idx = -1

        if idx >= 0:
            # Return the icon from the cache.
            iset = self._cache[idx]
        else:
            # Follow uic's naming convention.
            name = 'icon'
            idx = len(self._cache)

            if idx > 0:
                name += str(idx)

            icon = self._object_factory.createQObject("QIcon", name, (),
                    is_attribute=False)
            iset.set_icon(icon, self._qtgui_module)
            self._cache.append(iset)

        return iset.icon


class _IconSet(object):
    """An icon set, ie. the mode and state and the pixmap used for each."""

    def __init__(self, iconset, base_dir):
        """Initialise the icon set from an XML tag."""

        # Set the pre-Qt v4.4 fallback (ie. with no roles).
        self._fallback = self._file_name(iconset.text, base_dir)
        self._use_fallback = True

        # Parse the icon set.
        self._roles = {}

        for i in iconset:
            file_name = i.text
            if file_name is not None:
                file_name = self._file_name(file_name, base_dir)

            self._roles[i.tag] = file_name
            self._use_fallback = False

        # There is no real icon yet.
        self.icon = None

    @staticmethod
    def _file_name(fname, base_dir):
        """ Convert a relative filename if we have a base directory. """

        fname = fname.replace("\\", "\\\\")

        if base_dir != '' and fname[0] != ':' and not os.path.isabs(fname):
            fname = os.path.join(base_dir, fname)

        return fname

    def set_icon(self, icon, qtgui_module):
        """Save the icon and set its attributes."""

        if self._use_fallback:
            icon.addFile(self._fallback)
        else:
            for role, pixmap in self._roles.items():
                if role.endswith("off"):
                    mode = role[:-3]
                    state = qtgui_module.QIcon.Off
                elif role.endswith("on"):
                    mode = role[:-2]
                    state = qtgui_module.QIcon.On
                else:
                    continue

                mode = getattr(qtgui_module.QIcon, mode.title())

                if pixmap:
                    icon.addPixmap(qtgui_module.QPixmap(pixmap), mode, state)
                else:
                    icon.addPixmap(qtgui_module.QPixmap(), mode, state)

        self.icon = icon

    def __eq__(self, other):
        """Compare two icon sets for equality."""

        if not isinstance(other, type(self)):
            return NotImplemented

        if self._use_fallback:
            if other._use_fallback:
                return self._fallback == other._fallback

            return False

        if other._use_fallback:
            return False

        return self._roles == other._roles
