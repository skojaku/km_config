import os
from setuptools import setup
from setuptools.extension import Extension

this_dir = os.path.dirname(os.path.abspath(__file__))

#os.environ["CC"] = "g++"
#os.environ["CXX"] = "g++"

km_config = Extension(
    name='km_config',
    include_dirs=[os.path.join(this_dir, 'include')],
    sources=[
        os.path.join(this_dir, 'src', 'km_config.cpp')
    ],
    extra_compile_args=["-g0", "-O3", "-std=c++11", "-fopenmp"],
    language='c++11', 
)

setup(
    name='km_config',
    version='0.0.2',
    author='Sadamori Kojaku',
    description='Algorithm for finding multiple core-periphery pairs in networks',
    long_description='Algorithm for finding multiple core-periphery pairs in networks. Please cite  "S. Kojaku and N. Masuda, Core-periphery structure requires something else in the network" New J. Phys., 20, 043012 (2018)"',
    url='https://github.com/skojaku/km_config',
    classifiers=[
        # The list of PyPI classifiers
    ],
    ext_modules=[km_config],
    zip_safe=False,
    include_package_data=True,
)
