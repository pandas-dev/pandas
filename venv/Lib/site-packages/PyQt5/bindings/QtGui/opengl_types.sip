// This implements the typedefs for the OpenGL data types.
//
// Copyright (c) 2024 Riverbank Computing Limited <info@riverbankcomputing.com>
// 
// This file is part of PyQt5.
// 
// This file may be used under the terms of the GNU General Public License
// version 3.0 as published by the Free Software Foundation and appearing in
// the file LICENSE included in the packaging of this file.  Please review the
// following information to ensure the GNU General Public License version 3.0
// requirements will be met: http://www.gnu.org/copyleft/gpl.html.
// 
// If you do not wish to use this file under the terms of the GPL version 3.0
// then you may purchase a commercial license.  For more information contact
// info@riverbankcomputing.com.
// 
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


%If (PyQt_OpenGL)

typedef char GLchar;
typedef qint8 GLbyte;
typedef quint8 GLubyte;
typedef quint8 GLboolean;
typedef qint16 GLshort;
typedef quint16 GLushort;
typedef qint32 GLint;
typedef qint32 GLsizei;
typedef quint32 GLuint;
typedef quint32 GLenum;
typedef quint32 GLbitfield;
%If (PyQt_Desktop_OpenGL)
typedef quint64 GLuint64;       // This is in OpenGL ES v3.
typedef double GLdouble;
%End
typedef float GLfloat;
typedef float GLclampf;
typedef long GLintptr;
typedef long GLsizeiptr;

%End
