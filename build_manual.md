# Test Build
```bash
pip uninstall genmine
pip install -e .
```

# PyPI Build
In new conda environment:
```bash
conda create -n GenMine_build
conda activate GenMine_build
conda install pip
pip install hatch
hatch build
hatch publish  # Uses configured credentials
```

## Configure PyPI Token (Secure Method)
Instead of passing token as raw string, use hatch's secure configuration:
```bash
# Set username
hatch config set publish.index.pypi.username __token__

# Set token (will be stored securely in hatch config)
hatch config set publish.index.pypi.password your-pypi-token-here
```

Alternatively, use environment variables:
```bash
# Set environment variables
export HATCH_INDEX_USER=__token__
export HATCH_INDEX_AUTH=your-pypi-token-here

# Then publish without exposing token
hatch publish
```

Or create a `.pypirc` file in your home directory:
```ini
[pypi]
username = __token__
password = your-pypi-token-here
```
Then use: `hatch publish --repo pypi`

**Note**: The hatch config method is most secure as it stores credentials in hatch's encrypted config. Environment variables are good for CI/CD. The `.pypirc` file method matches the old twine workflow.
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
