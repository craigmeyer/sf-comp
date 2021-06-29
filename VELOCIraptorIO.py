"""
Routines for reading in VELOCIraptor hdf5 files

//todo: Extend to reading in data across multiple snapshots and split across multiple files
"""

import numpy as np
import pandas as pd
import h5py

# Constants for unit conversion
CM_TO_KPC = 3.085678e21
CMPERS_TO_KMPERS = 1e5
MSUN_TO_GRAM = 1.989e33

def readVELhalodata(base_fn, desiredFields = []):
    """
    Read in VELOCIraptor halo catalogue

    base_fn - base filename for VELOCIraptor output file
    desiredFields - desired fields to be read in, all fields read in if empty

    NOTE: Assumes VELOCIraptor halo catalogue is in hdf5 format
    """

    fn = base_fn + '.properties.0'

    halofile = h5py.File(fn, 'r')

    nhalo = np.uint64(halofile['Total_num_of_groups'])

    # Select the desired fields
    # If none are passed, read in all fields
    if len(desiredFields) > 0:
        fieldnames = desiredFields
    else:
        fieldnames = [str(fieldname) for fieldname in halofile.keys()]
        # Remove header info
        fieldnames.remove('File_id')
        fieldnames.remove('Num_of_files')
        fieldnames.remove('Num_of_groups')
        fieldnames.remove('Total_num_of_groups')
        fieldnames.remove('Configuration')
        fieldnames.remove('SimulationInfo')
        fieldnames.remove('UnitInfo')
    fieldtypes = [halofile[fieldname].dtype for fieldname in fieldnames]

    # Set up halo data dictionary
    halodata = {fieldnames[i]: np.zeros(nhalo, dtype = fieldtypes[i]) for i in range(len(fieldnames))}

    # Populate halo data dictionary
    for i in range(len(fieldnames)):
        ifield = fieldnames[i]
        halodata[ifield] = np.array(halofile[ifield])

    halofile.close()

    return halodata, nhalo

def readVELparticledata(base_fn):
    """
    Read in particle data from multiple VELOCIraptor output files

    base_fn - base filename for each VEL output file

    NOTE: Assumes the VELOCIraptor files are hdf5 format
    """

    gfn = base_fn + '.catalog_groups.0'
    pfn = base_fn + '.catalog_particles.0'
    upfn = base_fn + '.catalog_particles.unbound.0'

    gfile = h5py.File(gfn, 'r')
    nhalo = np.uint64(gfile['Num_of_groups'])
    n_inhalo = np.uint64(gfile['Group_Size']) # Total no. of parts in each halo
    offset = np.uint64(gfile['Offset'])
    uoffset = np.uint64(gfile['Offset_unbound'])
    gfile.close()

    pfile = h5py.File(pfn, 'r')
    pids = np.int64(pfile['Particle_IDs'])
    pfile.close()

    upfile = h5py.File(upfn, 'r')
    upids = np.int64(upfile['Particle_IDs'])
    unpart = len(upids) # Total number of unbound parts
    upfile.close()

    # Produce particle data dictionary
    partdata = {}
    partdata['Offset'] = offset
    partdata['Offset_unbound'] = uoffset
    partdata['Particle_IDs'] = [[] for ihalo in range(int(nhalo))]

    # Number of (bound/unbound) particles in each halo
    partdata['Npart'] = n_inhalo
    un_inhalo = np.zeros(nhalo, dtype = 'uint64') # No. of unbound parts in each halo
    for ihalo in range(int(nhalo - 1)):
        un_inhalo[ihalo] = uoffset[ihalo + 1] - uoffset[ihalo]
    un_inhalo[-1] = unpart - uoffset[-1]
    partdata['Npart_unbound'] = un_inhalo

    # PIDs for each halo - list all bound particles, then all unbound particles for each halo
    bn_inhalo = n_inhalo - un_inhalo # Bound particles in each halo
    for ihalo in range(int(nhalo)):
        partdata['Particle_IDs'][ihalo] = np.zeros(n_inhalo[ihalo], dtype = np.int64)
        partdata['Particle_IDs'][ihalo][:bn_inhalo[ihalo]] = pids[offset[ihalo]:offset[ihalo] + bn_inhalo[ihalo]]
        partdata['Particle_IDs'][ihalo][bn_inhalo[ihalo]:] = upids[uoffset[ihalo]:uoffset[ihalo] + un_inhalo[ihalo]]

    return partdata


def getUnitInfo(base_fn):
    """
    Extract unit information from simulation for future
    easy conversion to desired units

    Parameters:
    fn - Filename of VELOCIraptor unit file

    Returns:
    unitinfo - Array containing the unit information
    """

    fn = base_fn + '.units'

    unitfile = open(fn, 'r')

    # Construct unit information dictionary
    unitinfo = {}

    for line in unitfile:

        # Extract conversion value
        if line[0][0] == 'L':
            length_unit_to_kpc = np.float64(line.strip().split(' : ')[1].split(' # ')[0])
            unitinfo['UnitLength_in_cm'] = length_unit_to_kpc * CM_TO_KPC
        elif line[0][0] == 'V':
            vel_unit_to_km_per_s = np.float64(line.strip().split(' : ')[1].split(' # ')[0])
            unitinfo['UnitVelocity_in_cm_per_s'] = vel_unit_to_km_per_s * CMPERS_TO_KMPERS
        elif line[0][0] == 'M':
            mass_unit_to_solarmass = np.float64(line.strip().split(' : ')[1].split(' # ')[0])
            unitinfo['UnitMass_in_g'] = mass_unit_to_solarmass * MSUN_TO_GRAM
        else:
            continue
    unitfile.close()

    return unitinfo














