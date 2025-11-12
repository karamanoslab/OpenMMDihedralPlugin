from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
import subprocess, sysconfig, site, os, glob
import openmm

PLUGIN_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
PLUGIN_BUILD = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build"))

CONDA_PREFIX = os.environ.get("CONDA_PREFIX", "")
CONDA_LIB = os.path.join(CONDA_PREFIX, "lib")
CONDA_INCLUDE = os.path.join(CONDA_PREFIX, "include")
PYTHON_INCLUDE = sysconfig.get_paths()["include"]


class CustomBuildExt(build_ext):
    def run(self):
        swig_cmd = [
            "swig", "-python", "-c++",
            "-I./swig_lib",
            "./dihedralplugin/dihedralplugin.i"
        ]
        subprocess.check_call(swig_cmd)
        super().run()


class CustomBuildPy(build_py):
    def run(self):
        # Ensure SWIG-generated Python file exists before packaging
        swig_py = os.path.join("dihedralplugin", "dihedralplugin.py")
        if not os.path.exists(swig_py):
            subprocess.check_call(["swig", "-python", "-c++", "-I./swig_lib", "dihedralplugin/dihedralplugin.i"])
        super().run()


extension = Extension(
    name="dihedralplugin._dihedralplugin",
    sources=[
        "dihedralplugin/dihedralplugin_wrap.cxx",
        os.path.join(PLUGIN_SRC, "DihedralForce.cpp"),
        os.path.join(PLUGIN_SRC, "DihedralForceImpl.cpp"),
    ],
    
    swig_opts=[
        '-c++',
        '-python',
        '-I./swig_lib'
    ],
    
    include_dirs=[PYTHON_INCLUDE, CONDA_INCLUDE, PLUGIN_SRC],
    library_dirs=[CONDA_LIB, PLUGIN_BUILD ],
    libraries=["DihedralPlugin", "OpenMM"],
    extra_compile_args=["-fPIC"],
)

setup(
    name="dihedralplugin",
    version="1.0.0",
    author="Theodoros Karamanos",
    description="Custom CPU-only Dihedral Force plugin for OpenMM",
    packages=["dihedralplugin"],
    package_dir={"dihedralplugin":"dihedralplugin"}, 
    package_data={"dihedralplugin": ["*.so", "*.py"]},
    ext_modules=[extension],
    cmdclass={"build_ext": CustomBuildExt, "build_py": CustomBuildPy},
   # install_requires=["openmm>=8.0"],
    python_requires=">=3.7",
)