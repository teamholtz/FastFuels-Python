import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
  name = 'fastfuels',
  packages = ['fastfuels'],
  version = '0.5.5',
  license='GNU GPLv3',
  description = '3D fuelscapes for the contiguous US',
  long_description = README,
  long_description_content_type="text/markdown",
  author = 'Lucas Wells',
  author_email = 'lucas@holtzforestry.com',
  url = 'https://github.com/holtzforestry/FastFuels-Python',
  keywords = ['fire model', 'fuelscape', 'wildfire'],
  install_requires=[
          # pinning all dependencies to avoid version mismatch
          'colorcet==2.0.6',
          'fsspec==2021.11.0',
          # numcodecs no longer includes msgpack?
          'msgpack==1.0.2',
          'numpy==1.21.4',
          'pyvista==0.28.1',
          's3fs==2021.11.0',
          'scipy==1.7.2',
          'shapely>=1.7.1',
          'zarr==2.8.3'
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
