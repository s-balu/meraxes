# Gstar library dependencies
import os

deps = {}

deps['exec'] = {
    'mpicc' : '/usr/local/x86_64/gnu/openmpi-1.10.2-psm/bin/mpicc',
	'cc' : '/usr/local/intel-15.3.0/composer_xe_2015.3.187/bin/intel64/icc',
    'git' : '/usr/bin/git',
}

deps['gsl'] = {
    'inclp' :'/usr/local/x86_64/gnu/gsl-1.9/include',
    'libp'  :'/usr/local/x86_64/gnu/gsl-1.9/lib',
    'lib'   : ['gsl', 'gslcblas'],
}

deps['hdf5'] = {
    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/include',
    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.10.0-openmpi-1.10.2-psm/lib',
    'lib'   : ['hdf5', 'hdf5_hl', 'z'],
}

deps['gbpCode'] = {
    'inclp' :'/home/smutch/3rd_party/gbpCode/myInclude',
    'libp'  :'/home/smutch/3rd_party/gbpCode/myLib/mpi',
    'lib'   : 'gbpLib',
}

deps['fftw'] = {
    'inclp' :'/home/smutch/3rd_party/fftw-3.3.3/include',
    'libp'  :'/home/smutch/3rd_party/fftw-3.3.3/lib',
    'lib'   : ['fftw3f_mpi', 'fftw3f'],
}

deps['cmocka'] = {
    'inclp' : '/home/smutch/3rd_party/cmocka/include',
    'libp'  : '/home/smutch/3rd_party/cmocka/lib',
    'lib'   : 'cmocka',
}

deps['mlog'] = {
    'inclp' : os.path.join(os.path.dirname(os.path.realpath(__file__)), '../mlog'),
    'libp' : os.path.join(os.path.dirname(os.path.realpath(__file__)), '../mlog'),
    'lib'   : ['mlog'],
}
