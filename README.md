# FastFuels

3D fuelscapes for the contiguous US

## Quickstart

Import the module

```python
import fastfuels
```

Now open a `.fio` resource. If you have one locally, specify the path and file name.

```python
fio = fastfuels.open('../../FuelsIO/data/test/demo.fio')
```

Otherwise, you can connect to the cloud hosted demo.

```python
fio = fastfuels.open('remote')
```
