from setuptools import setup

REQUIRED_PACKAGES = [
    'gdal >= 2.2.4',
    'scipy >= 1.1.0',
    'matplotlib >= 2.2.2',
    'numpy >= 1.13.3',
    'google-api-python-client >= 1.7.3',
    'oauth2client >= 0.2.2'
]

setup(
    name='gee_tools',
    packages=['gee_tools'],
    description='A collection of utility functions relating to Google Earth Engine.',
    version='0.0',
    url='https://github.com/george-azzari/gee_tools/',
    author='George Azzari',
    author_email='',
    keywords=['GEE', 'Earth Engine', 'Google Earth Engine'],
    install_requires=REQUIRED_PACKAGES,
)
