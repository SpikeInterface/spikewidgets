import setuptools

d = {}
exec(open("spikewidgets/version.py").read(), None, d)
version = d['version']
pkg_name = "spikewidgets"

setuptools.setup(
    name=pkg_name,
    version=version,
    author="Jeremy Magland, Alessio Buccino",
    author_email="jmagland@flatironinstitute.org",
    description="Python widgets for use with spike sorting and electrophysiology",
    url="https://github.com/magland/spikewidgets",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'spiketoolkit',
        'matplotlib'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )
)
