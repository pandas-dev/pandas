#!/usr/bin/env python
from __future__ import division, print_function

fib3_f = '''
C FILE: FIB3.F
      SUBROUTINE FIB(A,N)
C
C     CALCULATE FIRST N FIBONACCI NUMBERS
C
      INTEGER N
      REAL*8 A(N)
Cf2py intent(in) n
Cf2py intent(out) a
Cf2py depend(n) a
      DO I=1,N
         IF (I.EQ.1) THEN
            A(I) = 0.0D0
         ELSEIF (I.EQ.2) THEN
            A(I) = 1.0D0
         ELSE
            A(I) = A(I-1) + A(I-2)
         ENDIF
      ENDDO
      END
C END FILE FIB3.F
'''

def source_func(ext, build_dir):
    import os
    from distutils.dep_util import newer
    target = os.path.join(build_dir, 'fib3.f')
    if newer(__file__, target):
        f = open(target, 'w')
        f.write(fib3_f)
        f.close()
    return [target]

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('gen_ext', parent_package, top_path)
    config.add_extension('fib3',
                         [source_func]
                         )
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
