import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
  name = 'fastfuels',
  packages = ['fastfuels'],
  version = '0.5.3',
  license='GNU GPLv3',
  description = '3D fuelscapes for the contiguous US',
  long_description = README,
  long_description_content_type="text/markdown",
  author = 'Lucas Wells',
  author_email = 'lucas@holtzforestry.com',
  url = 'https://github.com/holtzforestry/FastFuels-Python',
  keywords = ['fire model', 'fuelscape', 'wildfire'],
  install_requires=[
          'colorcet',
          'numpy',
          'pyvista',
          's3fs',
          'scipy',
          'shapely',
          'zarr',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
  scripts = ['fastfuels/fastfuels_create_index.py']
)
