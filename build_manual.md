# test build
```
pip uninstall GenMine
python setup.py clean --all
python setup.py install
```
# pypi build
In new conda environment
```
conda create -n GenMine_build
conda install pip
pip install twine
python setup.py bdist_wheel --universal
python setup.py sdist
twine upload dist/GenMine-{YOUR_VERSION}* 	// use current build number, don't forget *
```
# conda build
In new conda environment
```
conda create -n GenMine_condabuild
conda activate GenMine_condabuild
conda install -c conda-forge grayskull packaging -y
conda install conda-build conda-verify anaconda-client git -y
grayskull pypi GenMine
anaconda login
conda config --set anaconda_upload no
conda-build ./genmine -c conda-forge
anaconda upload {Build file location} // copy tar.bz2 file from the stdout of conda-build
```
