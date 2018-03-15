from distutils.core import setup, Extension
import numpy

# define the extension module
detection = Extension('detection', sources=['detection_py3_omp.c'],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args=['-O3','-fopenmp'],
                          extra_link_args=['-lgomp'])

p_detection = Extension('p_detection', sources=['p_detection_py3_omp.c'],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args=['-O3','-fopenmp'],
                          extra_link_args=['-lgomp'])


tt_processing = Extension('tt_processing', sources=['tt_processing_py3.c'],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args=['-O3','-fopenmp'],
                          extra_link_args=['-lgomp'])

LOC_STALTA = Extension('LOC_STALTA', sources=['stalta_py3.c'],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args=['-O3'])

DET_STALTA = Extension('DET_STALTA', sources=['stalta_det_py3.c'],
                          include_dirs=[numpy.get_include()],
                          extra_compile_args=['-O3'])


# run the setup
setup(ext_modules=[detection,p_detection,tt_processing, LOC_STALTA, DET_STALTA])

# python3 installation_detector.py build_ext --inplace
