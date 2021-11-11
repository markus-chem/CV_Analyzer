#!/usr/bin/env python
import json
import numpy as np
from setuptools import setup, find_packages, Extension


if __name__ == '__main__':
    """
    The first part compiles the broad package, the necessary
    information is given in the setup.json file.

    The second part compiles the files in the package written in Cython.
    As of now this is the cov.pyx file in the _helpers directory. Quite
    possibly you will need to change the c_flags.
    """

    with open('setup.json', 'r') as info:
        kwargs = json.load(info)
    setup(
        **kwargs
        )

    # assuming gcc, set (aggressive) optimization flags,
    # to be appended to the default compile line
    c_flags = []
    c_flags.append("-O3")
    c_flags.append("-ffast-math")
    c_flags.append("-fopenmp")
    c_flags.append("-march=native")
    c_flags.append("-fPIC")
    c_flags.append("-fopt-info")

    ld_flags = c_flags

    ext = Extension("cvpreparer._helpers.cov",
                    sources=["cvpreparer/_helpers/cov.pyx"],
                    extra_compile_args=c_flags,
                    extra_link_args=ld_flags,
                    include_dirs=[np.get_include()]
    )
    setup(
        **kwargs,
        ext_modules=[ext]
        )
