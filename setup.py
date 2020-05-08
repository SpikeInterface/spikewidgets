import setuptools

d = {}
exec(open("spikewidgets/version.py").read(), None, d)
version = d['version']
pkg_name = "spikewidgets"
long_description = open("README.md").read()

setuptools.setup(
    name=pkg_name,
    version=version,
    author="Jeremy Magland, Alessio Buccino",
    author_email="jmagland@flatironinstitute.org",
    description="Python widgets for spike sorting and electrophysiology",
    url="https://github.com/SpikeInterface/spikewidgets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'spiketoolkit',
        'spikecomparison',
        'matplotlib',
        'MEAutility>=1.4.6'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )
)
