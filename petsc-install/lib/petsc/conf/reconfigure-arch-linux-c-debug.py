#!/home/projects/ppc64le/python/2.7.12/bin/python
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    ')',
    '--download-fblaslapack',
    '--prefix=/ascldap/users/snikolo/move/DEMSI_newLAMMPS/petsc-install',
    'PETSC_ARCH=arch-linux-c-debug',
  ]
  configure.petsc_configure(configure_options)
