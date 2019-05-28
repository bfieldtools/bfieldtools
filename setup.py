from setuptools import setup, find_packages
import sys
import os
import numpy

#First, prevent accidental upload to PyPI
def forbid_publish():
    argv = sys.argv
    blacklist = ['register', 'upload']

    for command in blacklist:
        if command in argv:
            values = {'command': command}
            print('Command "%(command)s" has been blacklisted, exiting...' %
                    values)
            sys.exit(2)


def main():

    forbid_publish()

    setup(
        name="bfieldtools",
	description="Magnetic field modelling tools",
        long_description=codecs.open('README.md', encoding='utf8').read()
        version_format='{tag}.dev{commitcount}+{gitsha}',
        setup_requires=['setuptools-git-version'],
        install_requires=['numpy','scipy','matplotlib', 'mayavi', 'quadpy', 'trimesh', 'cvxopt', 'joblib', 'numba', 'multiprocessing'],
        packages=find_packages(),
	include_dirs=[numpy.get_include()],
        package_data={'bfieldtools': ['example_meshes/*.stl', 'example_meshes/*.obj']},
        classifiers=['Intended Audience :: Science/Research',
                     'Intended Audience :: Developers',
                     'Programming Language :: Python',
                     'Topic :: Software Development',
                     'Topic :: Scientific/Engineering'],
        platforms='any'
    )

### Authors, license etc. still need to be added


if __name__ == '__main__':
    main()