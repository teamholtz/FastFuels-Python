#!/usr/bin/env python3

"""
This script creates an index of FastFuels zarr files based on their extents.
"""

__author__     = "Daniel Crawl, UCSD"
__date__       = "29 July 2021"
__version__    = "0.5.3"
__maintainer__ = "Lucas Wells"
__email__      = "lucas@holtzforestry.com"
__status__     = "Prototype"

import argparse
import numcodecs
import os
import sys
import zarr

def parse_extent(extent, extent_fmt):
    if extent_fmt == '(x1, y1), (x2, y2)' and len(extent) == 4:
        return extent
    elif extent_fmt == '[[x1, y1], [x2, y2]]':
        return extent[0][0], extent[0][1], extent[1][0], extent[1][1]
    else:
        raise Exception(f'Unknown extent format: {extent_fmt}')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', nargs='+', required=True, help='Names of files to index.')
    parser.add_argument('-o', nargs=1, required=True, help='Name of index file to create.')
    parser.add_argument('-r', action='store_true', help='Use relative paths instead of absolute.')
    parser.add_argument('-v', action='store_true', help='Turn on verbosity')
    
    args = vars(parser.parse_args())
    
    fios = args['i']
    output = args['o'][0]
    relative = args['r']
    verbose = args['v']

    print(f'creating {output}')

    z = zarr.open(output, mode='w')
    group = z.create_group('index')

    fio_files = []
    for fio_name in sorted(fios):

        if fio_name[0] != '/' and not relative:
            print('WARNING: not absolute path {fio_name},')
    
        try:
            fio_file = zarr.open(fio_name)
    
            if relative:
                fio_name = os.path.relpath(fio_name, start=os.path.dirname(output))
    
            fio_files.append({
                'name': fio_name,
                'fio': fio_file
            })
        except:
            print(f'WARNING: Could not open {fio_name}')
            continue
    
    names = group.create_dataset('name', shape=(len(fio_files)), dtype=str)
    extents = group.create_dataset('extent', shape=(len(fio_files)), dtype=object, object_codec=numcodecs.MsgPack())
    extent_fmts = group.create_dataset('extent_fmt', shape=(len(fio_files)), dtype=str)
    
    min_x = 999999999
    min_y = 999999999
    max_x = -999999999
    max_y = -999999999
    
    i = 0
    for fio_entry in fio_files:
    
        fio_file = fio_entry['fio']
    
        if verbose:
            print(fio_entry['name'], fio_file.attrs['extent'])
    
        names[i] = fio_entry['name']
        extents[i] = fio_file.attrs['extent']
        extent_fmts[i] = fio_file.attrs['extent_fmt']
    
        x1, y1, x2, y2 = parse_extent(extents[i], extent_fmts[i])
    
        if x1 < min_x:
            min_x = x1
        if y1 < min_y:
            min_y = y1
        if x2 > max_x:
            max_x = x2
        if y2 > max_y:
            max_y = y2
    
    
        if i == 0:
            proj = fio_file.attrs['proj']
            res = fio_file.attrs['resolution']
            units = fio_file.attrs['units']
        elif proj != fio_file.attrs['proj']:
            raise Exception('Cannot handle different projections.')
        elif res != fio_file.attrs['resolution']:
            raise Exception('Cannot handle different resolutions.')
        elif units != fio_file.attrs['units']:
            Exception('Cannot handle different units.')
    
    
    
        i += 1
    
    #print(min_x,min_y,max_x,max_y)
    z.attrs['extent'] = [[min_x,min_y],[max_x,max_y]]
    z.attrs['extent_fmt'] = '[[x1, y1], [x2, y2]]'
    z.attrs['proj'] = proj
    z.attrs['resolution'] = res
    z.attrs['units'] = units
    z = None
