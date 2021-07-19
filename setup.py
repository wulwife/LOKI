import os
import numpy
from setuptools import setup, find_packages, Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    required_list = f.read().splitlines()

# ===================================  DEFINE COMPILER PATH !!!
# os.environ["CC"] = 'gcc-mp-7'
os.environ["CC"] = "gcc"


# ===================================  C - COMPILING
location = Extension('location',
                     sources=['loki/src_c/location_py3_omp.c'],
                     include_dirs=[numpy.get_include()],
                     extra_compile_args=['-O3', '-fopenmp'],
                     extra_link_args=['-lgomp'])

location_t0 = Extension('location_t0',
                     sources=['loki/src_c/location_t0_py3_omp.c'],
                     include_dirs=[numpy.get_include()],
                     extra_compile_args=['-O3', '-fopenmp'],
                     extra_link_args=['-lgomp'])

detection = Extension('detection',
                      sources=['loki/src_c/detection_py3_omp.c'],
                      include_dirs=[numpy.get_include()],
                      extra_compile_args=['-O3', '-fopenmp'],
                      extra_link_args=['-lgomp'])

tt_processing = Extension('tt_processing',
                          sources=['loki/src_c/tt_processing_py3.c'],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args=['-O3', '-fopenmp'],
                          extra_link_args=['-lgomp'])

LOC_STALTA = Extension('LOC_STALTA',
                       sources=['loki/src_c/stalta_py3.c'],
                       include_dirs=[numpy.get_include()],
                       extra_compile_args=['-O3'])

DET_STALTA = Extension('DET_STALTA',
                       sources=['loki/src_c/stalta_det_py3.c'],
                       include_dirs=[numpy.get_include()],
                       extra_compile_args=['-O3'])


# ===================================  SETUP
setup(
    name="loki",
    version="1.0.0",
    author="Francesco Grigoli",
    author_email="fsco.grigoli@gmail.com",
    description="Location of seismic events through traveltime staking",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wulwife/LOKI",
    python_requires='>=3.6',
    install_requires=required_list,
    packages=find_packages(),
    package_data={"loki": ['src_c/*.c']},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix, MacOS",
        "Intended Audience :: Science/Research",
    ],
    ext_modules=[location,
                 location_t0,
                 detection,
                 tt_processing,
                 LOC_STALTA,
                 DET_STALTA]
)


# =================================== HELP

# In this example, setup() is called with additional meta-information,
# which is recommended when distribution packages have to be built.
# For the extension itself, it specifies preprocessor defines,
# include directories, library directories, and libraries.
# Depending on the compiler, distutils passes this information in
# different ways to the compiler.
# For example, on Unix, this may result in the compilation commands

# ```
# gcc -DNDEBUG -g -O3 -Wall -Wstrict-prototypes -fPIC -DMAJOR_VERSION=1 \
#     -DMINOR_VERSION=0 -I/usr/local/include \
#     -I/usr/local/include/python2.2 -c demo.c \
#     -o build/temp.linux-i686-2.2/demo.o

# gcc -shared build/temp.linux-i686-2.2/demo.o -L/usr/local/lib -ltcl83
#     -o build/lib.linux-i686-2.2/demo.so
# ```

# These lines are for demonstration purposes only; distutils users
# should trust that distutils gets the invocations right.

# The manifest template has one command per line, where each command specifies a set of files to include or exclude from the source distribution. For an example, again we turn to the Distutils’ own manifest template:

# include *.txt
# recursive-include examples *.txt *.py
# prune examples/sample?/build

# The meanings should be fairly clear: include all files in the
# distribution root matching *.txt, all files anywhere under the examples
# directory matching *.txt or *.py, and exclude all directories matching
# examples/sample?/build. All of this is done after the standard include
# set, so you can exclude files from the standard set with explicit
# instructions in the manifest template. (Or, you can use the --no-defaults
# option to disable the standard set entirely.)
# There are several other commands available in the manifest template
# mini-language; see section Creating a source distribution:
# the sdist command.
