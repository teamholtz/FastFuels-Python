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

You can also view metadata for resolution and units

```python
>>> print(fio.res)
(1,1,1)
>>> print(fio.units)
'meters'
```

### Spatial queries

You can perform spatial queries by specifying geographic coordinates in decimal degrees
and a radius in meters. The radius parameter defines the size of the bounding square in which
fuels are queried.

```python
# this command will return a square kilometer of fuels (radius=500 meters)
roi = fio.query(-122.191, 41.208, 500)
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
