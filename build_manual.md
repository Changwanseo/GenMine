	# test build
python setup.py install

# pypi build
python setup.py bdist_wheel
twine upload dist/GenMine-*.*.*.*-py3-none-any.whl 	// use current build number

# conda build
conda-build ./conda/ -c bioconda
