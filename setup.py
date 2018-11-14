from setuptools import setup, find_packages

REQUIRED_PACKAGES = [
    'gdal >= 1.10.1',
    'scipy >= 1.1.0',
    'matplotlib >= 2.2.2',
    'numpy >= 1.13.3',
    'google-api-python-client >= 1.7.3',
    'oauth2client >= 0.2.2'
]

setup(
    name='gee_tools',
    packages=find_packages(exclude=['*.csv', '*.ipynb']),
    description='A collection of utility functions relating to Google Earth Engine.',
    version='0.0.3',
    url='https://github.com/AtlasAIPBC/gee_tools.git',
    author='George Azzari',
    author_email='',
    keywords=['GEE', 'Earth Engine', 'Google Earth Engine'],
    install_requires=REQUIRED_PACKAGES,
    python_requires="~=2.7"
)
