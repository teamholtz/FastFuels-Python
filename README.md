# FastFuels

3D fuelscapes for the contiguous US

## Install

You can install FastFuels through the Python Package Index.

```
pip install fastfuels
```

## Quickstart

### Connecting to a `.fio` resource

A `.fio` resource is a directory-in-file object where important metadata and fuel arrays are stored. Start by importing the FastFuels module and open a `.fio` resource. If you have one locally, specify the path and file name.

```python
>>> import fastfuels
>>> fio = fastfuels.open('./demo.fio')
```

Otherwise, you can connect to the cloud hosted demo.

```python
>>> fio = fastfuels.open('remote')
connecting to remote FIO server...
```

### Explore the metadata

Let's take a look at some metadata. You can get the extent of the data in geographic coordinates (longitude and latitude) or in projected coordinates by changing the `mode` argument.

```python
>>> print(fio.get_extent(mode='geographic'))
(-120.73665218037868, 38.93933418427242, -120.6511979123941, 38.90135366961076)
>>> print(fio.get_extent(mode='projected'))
(-2100315.0, 2043015.0, -2094315.0, 2037015.0)
```

And the projection system is stored in the `proj` attribute.

```python
>>> print(fio.proj)
```

You can also view metadata for resolution, dimensions and units

```python
>>> print(fio.res)
(1,1,1)
>>> print(fio.units)
'meters'
>>> print(fio.dim)
(6000, 6000, 100)
>>> print(fio.dim_fmt)
'x,y,z'
```

### Spatial queries

You can perform spatial queries by supplying the coordinates of a bounding box in three different modes. Relative queries reference the datasets relative to the upper left-hand corner of the horizontal dimensions. The following will extract a 300x300 meter region of interest (ROI).

```python
roi = fio.query((2000, 2000), (2300, 2300), mode='relative')
```

Queries can also be performed using geographic or projected coordinates to define the bounding box.

```python
roi = fuels.query((-2098000, 2039000), (-2097700, 2038700), mode='projected')
roi = fuels.query((-120.71, 38.93), (-120.705, 38.9275), mode='geographic')
```

### Viewing fuels in 3D

Fuel parameter arrays can be viewed interactively in 3D. To see the available parameters run

```python
print(roi.get_properties())
```

Then specify one of the properties in the `view()` method on the `roi` object.

```python
roi.view('sav')
```

![FastFuels SAV](https://storage.googleapis.com/public-assests/fastfuels_sav.png)

### Writing fire model input files

With the `roi` object, you can write input files for various fire models. Here,
you may also decrease the resolution to save computation.

```python
roi.write('./outputs', model='quicfire', res_xyz=[2,2,1])
```
