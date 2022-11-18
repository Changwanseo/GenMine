# test build
```
pip uninstall GenMine
python setup.py clean --all
python setup.py install
```
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
conda install -c conda-forge grayskull
conda install -c conda-forge packaging
conda install conda-build
conda install conda-verify
conda install anaconda-client
conda install git
grayskull pypi GenMine
anaconda login
conda config --set anaconda_upload no
conda-build ./genmine -c conda-forge
anaconda upload {Build file location}
```
