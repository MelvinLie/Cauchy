import os
import numpy
import platform
import shutil
from distutils.core import Distribution, Extension

from Cython.Build import build_ext, cythonize

# ========================================
# Build file written by M. Liebsch
# 
# we wrap the cbem code using cython
# ========================================



# compiler arguments are system dependent
if platform.system() == 'Windows':
    compile_args = ['/std:c++17', '/O2', '-openmp']
    extra_link_args = []
elif platform.system() == 'Linux':
    compile_args = ['-std=c++17', '-O3', '-fopenmp']
    extra_link_args = ['-fopenmp']

cython_dir = "cauchy_c"
extension = Extension(
    "cauchy_c",
    [
        os.path.join(cython_dir, "interface.pyx"),
    ],
    extra_compile_args=compile_args,
    include_dirs=[numpy.get_include()],
    extra_link_args=extra_link_args,
)

ext_modules = cythonize([extension], include_path=[cython_dir,
                                                   os.path.join("cauchy_c", "c-algorithms")])
dist = Distribution({"ext_modules": ext_modules})
cmd = build_ext(dist)
cmd.ensure_finalized()
cmd.run()

for output in cmd.get_outputs():
    relative_extension = os.path.relpath(output, cmd.build_lib)
    shutil.copyfile(output, os.path.join("cauchy", relative_extension))