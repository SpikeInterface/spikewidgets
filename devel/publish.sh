rm -rf build/ dist/ ml_ephys.egg-info/
python3 setup.py sdist bdist_wheel
twine upload dist/*
