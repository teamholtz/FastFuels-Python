# FastFuels

3D fuelscapes for the contiguous US

## Quickstart

### Connecting to a `.fio` resource

A `.fio` resource is a directory-in-file object where important metadata and fuel arrays are stores. Start by importing the FastFuels module and open a `.fio` resource. If you have one locally, specify the path and file name.

```python
>>> import fastfuels
>>> fio = fastfuels.open('../../FuelsIO/data/test/demo.fio')
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
PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG
","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETE
R["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]
```
