import setuptools

d = {}
exec(open("spikewidgets/version.py").read(), None, d)
version = d['version']
pkg_name = "spikewidgets"
long_description = open("README.md").read()

setuptools.setup(
    name=pkg_name,
    version=version,
    author="Alessio Buccino, Cole Hurwitz, Samuel Garcia, Jeremy Magland, Matthias Hennig",
    author_email="jmagland@flatironinstitute.org",
    description="Python widgets for spike sorting and electrophysiology",
    url="https://github.com/SpikeInterface/spikewidgets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'spikeextractors>=0.9.7',
        'spiketoolkit>=0.7.6',
        'spikecomparison>=0.3.3',
        'matplotlib',
        'MEAutility>=1.4.8'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )
)
