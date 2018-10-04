import setuptools

pkg_name="spikewidgets"

import unittest
def my_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite

setuptools.setup(
    name=pkg_name,
    version="0.1.7",
    author="Jeremy Magland",
    author_email="jmagland@flatironinstitute.org",
    description="Python widgets for use with spike sorting and electrophysiology",
    url="https://github.com/magland/spikewidgets",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'ipython',
        'vdom',
        'ipywidgets',
        'matplotlib',
        'spikeinterface'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ),
    test_suite='setup.my_test_suite'
)
