"""
Routines for reading in VELOCIraptor hdf5 files

//todo: Extend to reading in data across multiple snapshots and split across multiple files
"""

import numpy as np
import pandas as pd
import h5py
import velociraptor_python_tools as vpt

# Constants for unit conversion
# CM_TO_MPC = 3.085678e24
# CMPERS_TO_KMPERS = 1e5
# MSUN_TO_GRAM = 1.989e33

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
        filednames.remove('Num_of_groups')
        fieldnames.remove('Total_num_of_groups')
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
    partdata['Npart'] = np.zeros(nhalo, dtype = np.uint64)
    partdata['Npart_unbound'] = np.zeros(nhalo, dtype = np.uint64)
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














