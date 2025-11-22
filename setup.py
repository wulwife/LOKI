import os
import sys
import platform
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

# ===== Long description / requirements =====
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
    required_list = f.read().splitlines()

# ====== Compilatore (opzionale) ======
# os.environ["CC"] = "gcc"  # lasciare commentato se non necessario

# ====== OpenMP flags cross-platform ======
def _omp_flags():
    """
    Restituisce (extra_compile_args, extra_link_args) per OpenMP
    a seconda della piattaforma/compilatore.
    """
    sysname = platform.system()
    if sysname == "Darwin":
        # Clang su macOS richiede libomp (brew install libomp)
        # e flag diversi
        return (["-O3", "-Xpreprocessor", "-fopenmp"],
                ["-lomp"])
    elif sysname == "Windows":
        # MSVC
        return (["/O2", "/openmp"], [])
    else:
        # Linux/Unix con gcc
        return (["-O3", "-fopenmp"],
                ["-lgomp"])

omp_compile_args, omp_link_args = _omp_flags()

# ====== Estensioni C (senza numpy.get_include qui!) ======
location = Extension(
    "location",
    sources=["loki/src_c/location_py3_omp.c"],
    include_dirs=[],
    extra_compile_args=omp_compile_args,
    extra_link_args=omp_link_args,
)

location_t0 = Extension(
    "location_t0",
    sources=["loki/src_c/location_t0_py3_omp.c"],
    include_dirs=[],
    extra_compile_args=omp_compile_args,
    extra_link_args=omp_link_args,
)

voronoi_loc = Extension(
    "voronoi_loc",
    sources=["loki/src_c/voronoi_loc.c"],
    include_dirs=[],
    extra_compile_args=omp_compile_args,
    extra_link_args=omp_link_args,
)

detection = Extension(
    "detection",
    sources=["loki/src_c/detection_py3_omp.c"],
    include_dirs=[],
    extra_compile_args=omp_compile_args,
    extra_link_args=omp_link_args,
)

tt_processing = Extension(
    "tt_processing",
    sources=["loki/src_c/tt_processing_py3.c"],
    include_dirs=[],
    extra_compile_args=omp_compile_args,
    extra_link_args=omp_link_args,
)

LOC_STALTA = Extension(
    "LOC_STALTA",
    sources=["loki/src_c/stalta_py3.c"],
    include_dirs=[],
    extra_compile_args=["-O3"],
    extra_link_args=[],
)

DET_STALTA = Extension(
    "DET_STALTA",
    sources=["loki/src_c/stalta_det_py3.c"],
    include_dirs=[],
    extra_compile_args=["-O3"],
    extra_link_args=[],
)

ext_modules = [
    location,
    location_t0,
    voronoi_loc,
    detection,
    tt_processing,
    LOC_STALTA,
    DET_STALTA,
]

# ====== build_ext custom: aggiunge numpy.get_include() al volo ======
class BuildExtWithNumpy(build_ext):
    def finalize_options(self):
        super().finalize_options()
        # importa numpy SOLO ora (dopo che pip l'ha installato via pyproject)
        import numpy
        np_inc = numpy.get_include()
        for ext in self.extensions:
            if np_inc not in ext.include_dirs:
                ext.include_dirs.append(np_inc)

# ====== setup() ======
setup(
    name="loki",
    version="1.0.0",
    author="Francesco Grigoli",
    author_email="fsco.grigoli@gmail.com",
    description="Location of seismic events through traveltime stacking",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wulwife/LOKI",
    python_requires=">=3.8",
    install_requires=required_list,
    packages=find_packages(),
    package_data={"loki": ["src_c/*.c"]},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix, MacOS",
        "Intended Audience :: Science/Research",
    ],
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExtWithNumpy},
)
