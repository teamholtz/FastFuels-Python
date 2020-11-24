from distutils.core import setup

setup(
  name = 'FastFuels',
  packages = ['fastfuels'],
  version = '0.1',
  license='GNU GPLv3',
  description = '3D fuelscapes for the contiguous US'
  author = 'Lucas Wells',
  author_email = 'lucas@holtzforestry.com',
  url = 'https://github.com/holtzforestry/FastFuels',
  download_url = 'https://github.com/holtzforestry/FastFuels/archive/v0.1-alpha.tar.gz',
  keywords = ['fire model', 'fuelscape', 'wildfire'],
  install_requires=[
          'colorcet',
          'gcsfs',
          'numpy',
          'pyvista',
          'scipy',
          'zarr',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GNU GPLv3 License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
