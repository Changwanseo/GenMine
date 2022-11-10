# test build
python setup.py install

# pypi build

In new conda environment
```
conda install pip
pip install twine
python setup.py bdist_wheel --universal
python setup.py sdist
twine upload dist/GenMine-{YOUR_VERSION}* 	// use current build number
```
# conda build

In new conda environment
```
conda install conda-build
conda install conda-verify
conda install git
conda skeleton pypi GenMine --version {YOUR_VERSION}
#conda-build ./conda/ -c bioconda
```
