"""
Routines for reading in SubFind and Gadget4 hdf5 files

//todo: Extend to reading in data across multiple snapshots and split across multiple files
"""

import numpy as np
import h5py

def readSnapshot(opt, fn, desiredFields):
    ''' 
    Read in Gadget4 snapshot

    opt - Configuration options
    fn - snapshot filename
    desiredFields - any desired fields to be read in 
    '''

    snapfile = h5py.File(fn, 'r')

    # Read in desired fields from the SubFind halo catalogue
    snapdata = {}

    for gname, dsnames in desiredFields.items():
        group = snapfile[gname]
        for dsname in dsnames:
            full_name = '%s/%s' % (gname, dsname)
            snapdata[full_name] = np.array(group[dsname])

    return list(snapdata.values())


def loadSFhaloCatalogue(fn, desiredFields = []):
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

    if len(desiredFields) > 0:
        for gname, dsnames in desiredFields.items():
            group = halofile[gname]
            for dsname in dsnames:
                full_name = '%s/%s' % (gname, dsname)
                halodata[full_name] = group[dsname][()]
    else:
        groups = list(halofile.keys())
        groups.remove('Config')
        groups.remove('Header')
        groups.remove('IDs')
        groups.remove('Parameters')
        for gname in groups:
            group = halofile[gname]
            dsnames = list(halofile[gname].keys())
            for dsname in dsnames:
                full_name = '%s/%s' % (gname, dsname)
                halodata[full_name] = group[dsname][()]

    halofile.close()

    return halodata, nhalo, nsub

def getUnitInfo(fn):
    """
    Extract unit information from simulation for future
    easy conversion to desired units

    Parameters:
    fn - Filename of Gadget4 snapshot containing unit information

    Returns:
    unitinfo - Array containing the unit information
    """

    snapfile = h5py.File(fn, 'r')

    # Read in unit information from parameters
    params = snapfile['Parameters'].attrs
    length_unit = params['UnitLength_in_cm']
    mass_unit = params['UnitMass_in_g']
    vel_unit = params['UnitVelocity_in_cm_per_s']

    snapfile.close()

    # Construct unit information dictionary
    unitinfo = {}
    unitinfo['UnitLength_in_cm'] = length_unit
    unitinfo['UnitMass_in_g'] = mass_unit
    unitinfo['UnitVelocity_in_cm_per_s'] = vel_unit

    return unitinfo

















