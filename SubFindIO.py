"""
Routines for reading in SubFind and Gadget4 hdf5 files

//todo: Extend to reading in data across multiple snapshots and split across multiple files
"""

import numpy as np
import h5py

# Constants for unit conversion
# CM_TO_MPC = 3.085678e24
# CMPERS_TO_KMPERS = 1e5
# MSUN_TO_GRAM = 1.989e33

def readSnapshot(opt, fn, desiredFields):
    ''' 
    Read in Gadget4 snapshot

    opt - Configuration options
    fn - snapshot filename
    desiredFields - any desired fields to be read in 
    '''

    snapfile = h5py.File(fn, 'r')

    # Read in parameters for unit conversion
    # params = snapfile['Parameters'].attrs
    # length_unit = params['UnitLength_in_cm']
    # vel_unit = params['UnitVelocity_in_cm_per_s']

    # Read in desired fields from the SubFind halo catalogue
    snapdata = {}

    for gname, dsnames in desiredFields.items():
        group = snapfile[gname]
        for dsname in dsnames:
            full_name = '%s/%s' % (gname, dsname)
            snapdata[full_name] = np.array(group[dsname])

    return list(snapdata.values())


def loadSFhaloCatalogue(fn, desiredFields):
    ''' 
    Load in SubFind halo catalogue

    fn - filename of halo catalogue
    desiredFields - any desired fields to be read in
    '''

    halofile = h5py.File(fn, 'r')

    # Read in header attributes
    header = halofile['Header'].attrs
    nhalo = header['Ngroups_ThisFile']
    nsub = header['Nsubhalos_ThisFile']

    # Read in desired fields from the SubFind halo catalogue
    halodata = {}
    for gname, dsnames in desiredFields.items():
        group = halofile[gname]
        for dsname in dsnames:
            full_name = '%s/%s' % (gname, dsname)
            halodata[full_name] = group[dsname][()]

    halofile.close()

    return list(halodata.values()), nhalo, nsub















