"""
Main file for crossmatching two (sub)halo catalogues from separate structure finder codes.

Given two halo catalogues run on the same simulation snapshot and a configuration file, 
this script will find all matching structures in the two catalogues, save separate 
crossmatched catalogues for halos and subhalos, and print some simple summary statistics.

The final catalogues lists the ID of each (sub)halo from the reference catalogue and the ID 
of its corresponding match in the target catalogue.

//todo: extend to reading in catalogues across multiple snapshots or catalogues split across multiple files
"""

import sys
import h5py
import os
from ReadConfig import  CrossmatchOptions
import SubFindIO
import VELOCIraptorIO
import pandas as pd
import numpy as np

##################################################

def makeSFtables(opt, labelsDict):
    ''' 
    Generate a pandas table from SubFind halo catalogue

    opt - Configuration options
    labelsDict - dictionary specifying the column names in the table

    //todo: Check for file existence when saving to csv
    '''

    snap_fn = opt.snapdir + '/snapshot_%03d.hdf5' % opt.snap
    fields = {'PartType1': ('ParticleIDs',)}
    pids = SubFindIO.readSnapshot(opt, snap_fn, fields)[0]

    halo_fn = opt.SFdir + '/fof_subhalo_tab_%03d.hdf5' % opt.snap
    # Relevant fields to be read in 
    fields = {'Group': ('GroupLen', 'GroupOffsetType', 'GroupFirstSub'),
              'Subhalo': ('SubhaloLen', 'SubhaloOffsetType')}

    # Read in the relevant fields from the SubFind halo catalogue
    halodata, nhalo, nsub = SubFindIO.loadSFhaloCatalogue(halo_fn, fields)
    npart, offset, firstSub, npart_sub, offset_sub = halodata.values()
    offset = offset[:,1]
    offset_sub = offset_sub[:,1]

    # Create contiguous array of host group for each particle and PIDs for each halo
    hostHalo = np.empty(nhalo, dtype = 'object')
    haloPIDs = np.empty(nhalo, dtype = 'object')
    linds = offset
    uinds = offset + npart
    for ihalo in range(nhalo):
        hostHalo[ihalo] = np.array([ihalo] * npart[ihalo])
        haloPIDs[ihalo] = np.array([part for part in pids[linds[ihalo]:uinds[ihalo]]])
    hostHalo = np.concatenate(hostHalo)

    # Create contiguous array of host subhalo for each particle
    hostSubhalo = np.empty(nsub, dtype = 'object')
    for isub in range(int(nsub - 1)):
        # -1 entries are for background particles that are in the FOF group but not found in substructure
        hostSubhalo[isub] = np.array(([isub] * npart_sub[isub]) + ([-1] * ((offset_sub[isub + 1] - offset_sub[isub]) - npart_sub[isub])))
    # Add in the particles in the final subhalo
    hostSubhalo[int(nsub - 1)] = np.array([int(nsub - 1)] * (offset_sub[-1] - offset_sub[-2]))
    hostSubhalo = np.concatenate(hostSubhalo)

    haloPIDs = np.concatenate(haloPIDs)

    # Construct the pandas data tables - one for all structures, one for just subhaloes
    cols = {labelsDict['PID']: haloPIDs, labelsDict['HostHalo']: hostHalo, labelsDict['HostSubhalo']: hostSubhalo}
    table = pd.DataFrame(data = cols)
    table_sub = table.loc[~np.in1d(np.array(table[labelsDict['HostSubhalo']]), firstSub)]

    # Save the data products as csv files
    table.to_csv(opt.SFoutfilename + '_halo.csv', index = False)
    table_sub.to_csv(opt.SFoutfilename + '_sub.csv', index = False)

    return table, table_sub

def makeVELtables(opt, labelsDict):
    ''' 
    Generate a pandas table from VELOCIraptor halo catalogue

    opt - Configuration options
    labelsDict - dictionary specifying the column names in the table

    //todo: Check for file existence when saving to csv
    '''

    base_fn = opt.VELdir + '/snapshot_%03d.VELOCIRAPTOR' % opt.snap

    # Read in the VELOCIraptor halo and particle properties
    halodata, _ = VELOCIraptorIO.readVELhalodata(base_fn, desiredFields = ['hostHaloID'])
    partdata = VELOCIraptorIO.readVELparticledata(base_fn)

    # Select the relevant fields
    pids = partdata['Particle_IDs']
    npart = partdata['Npart']
    hostID = halodata['hostHaloID']
    nhalo = np.where(hostID == -1)[0].size
    nsub = hostID.size - nhalo

    # Create contiguous array of host group for each particle, handling field haloes and subhaloes separately
    hostHalo = np.empty(nhalo + nsub, dtype = 'object')
    for ihalo in range(nhalo):
        hostHalo[ihalo] = np.array([ihalo] * npart[ihalo])
    for isub in range(nhalo, nhalo + nsub):
        hostHalo[isub] = np.array([hostID[isub] - 1] * npart[isub])
    hostHalo = np.concatenate(hostHalo)

    # Create contiguous array of host subhalo for each particle, handling field haloes and subhaloes separately
    hostSubhalo = np.empty(nhalo + nsub, dtype = 'object')
    for ihalo in range(nhalo):
        # Particles in the (background) field haloes register no subhalo ID (-1)
        hostSubhalo[ihalo] = np.array([-1] * npart[ihalo])
    for isub in range(nhalo, nhalo + nsub):
        # Subhalo IDs start from 0 at the halo after the final field halo
        hostSubhalo[isub] = np.array([isub - nhalo] * npart[isub])
    hostSubhalo = np.concatenate(hostSubhalo)

    pids = np.concatenate(pids)

    # Construct the pandas data table - one for all structures, one for just subhaloes
    cols = {labelsDict['PID']: pids, labelsDict['HostHalo']: hostHalo, labelsDict['HostSubhalo']: hostSubhalo}
    table = pd.DataFrame(data = cols)
    table_sub = table.loc[table[labelsDict['HostSubhalo']] > -1]

    # Save the data products as csv files
    table.to_csv(opt.VELoutfilename + '_halo.csv', index = False)
    table_sub.to_csv(opt.VELoutfilename + '_sub.csv', index = False)

    return table, table_sub

def crossmatchCatalogues(opt, ref_table, target_table, ref_subtable, target_subtable, ref_finder, target_finder):
    '''
    Crossmatch two halo catalogues to find common members in both
    opt - Configuration options
    ref_table - PIDs, host halo for each particle in reference catalogue
    target_table - same as ref_table, but for target catalogue
    ref_subtable - PIDs, host subhalo for each particle in reference catalogue
    target-subtable - same as ref_subtable, but for target catalogue
    ref_finder - Name of reference structure finder
    target_finder - Name of target structure finder

    //todo: Check for file existence when saving to csv
    '''

    # Join the two catalogues, matching rows based on particle ID
    joinTable = pd.merge(ref_table, target_table, on = 'PID')

    # For each reference catalogue FOF halo, match it with the target catalogue FOF halo which 
    # contains the most common particles
    matchTable = joinTable.groupby(ref_table.columns[1])[target_table.columns[1]].apply(
        lambda x: x.value_counts().index[0]).reset_index()
    numMatches = len(matchTable)

    print('There are %s FOF haloes matched between the %s and %s catalogues' % (numMatches, ref_finder, target_finder))
    nhalo_ref = pd.unique(ref_table[ref_table.columns[1]]).size
    nhalo_target = pd.unique(target_table[target_table.columns[1]]).size
    print('These catalogues contained %s and %s FOF haloes, respectively' % (nhalo_ref, nhalo_target))

    # Repeat the matching scheme for subhaloes
    # We first only include rows that have non-negative subhalo IDs in both catalogues
    # This ensures we are only matching across particles in subhaloes
    joinTable_sub = pd.merge(ref_subtable, target_subtable)
    matchTable_sub = joinTable_sub.groupby(ref_subtable.columns[2])[target_subtable.columns[2]].apply(
        lambda x: x.value_counts().index[0]).reset_index()
    numMatches_sub = len(matchTable_sub)

    print('There are %s subhaloes matched between the %s and %s catalogues' % (numMatches_sub, ref_finder, target_finder))
    nsubhalo_ref = pd.unique(ref_subtable[ref_subtable.columns[2]]).size
    nsubhalo_target = pd.unique(target_subtable[target_subtable.columns[2]]).size
    print('These catalogues contained %s and %s subhaloes, respectively' % (nsubhalo_ref, nsubhalo_target))

    # Save the crossmatched catalogues as csv files
    matchTable.to_csv(opt.catoutfilename + '_halo.csv', index = False)
    matchTable_sub.to_csv(opt.catoutfilename + '_sub.csv', index = False)


def main(ref_finder, target_finder, config_file):

    # SubFind and VELOCIraptor - may be expanded in the future (most likely with HBT)
    available_finders = ['SF-HBT', 'VEL'] 

    # Check valid structure finders have been specified
    if ref_finder not in available_finders:
        raise SystemExit('Structure finder %s not recognised.' % ref_finder, 'Please choose from ', available_finders)
    if target_finder not in available_finders:
        raise SystemExit('Structure finder %s not recognised.' % target_finder, 'Please choose from ', available_finders)

    opt = CrossmatchOptions(ref_finder, target_finder, config_file)

    # Set up notation for the two finders (used for column names in tables)
    if ref_finder == 'SF-HBT':

        print('Reading in SubFind halo catalogue as reference catalogue')

        # Use the SubFind terminology
        labelsDict_ref = {}
        labelsDict_ref['PID'] = 'PID'
        labelsDict_ref['HostHalo'] = 'SF-HBT FOF Group'
        labelsDict_ref['HostSubhalo'] = 'SF-HBT Subhalo'

        ref_table, ref_subtable = makeSFtables(opt, labelsDict_ref)
    elif ref_finder == 'VEL':

        print('Reading in VELOCIraptor halo catalogue as reference catalogue')

        # Use the VELOCIraptor terminology
        labelsDict_ref = {}
        labelsDict_ref['PID'] = 'PID'
        labelsDict_ref['HostHalo'] = 'VEL Field Halo'
        labelsDict_ref['HostSubhalo'] = 'VEL Subhalo'

        ref_table, ref_subtable = makeVELtables(opt, labelsDict_ref)

    if target_finder == 'SF-HBT':

        print('Reading in SubFind halo catalogue as target catalogue')

        labelsDict_target = {}
        labelsDict_target['PID'] = 'PID'
        labelsDict_target['HostHalo'] = 'SF-HBT FOF Group'
        labelsDict_target['HostSubhalo'] = 'SF-HBT Subhalo'

        target_table, target_subtable = makeSFtables(opt, labelsDict_target)
    elif target_finder == 'VEL':

        print('Reading in VELOCIraptor halo catalogue as target catalogue')

        labelsDict_target = {}
        labelsDict_target['PID'] = 'PID'
        labelsDict_target['HostHalo'] = 'VEL Field Halo'
        labelsDict_target['HostSubhalo'] = 'VEL Subhalo'

        target_table, target_subtable = makeVELtables(opt, labelsDict_target)

    # Perform the crossmatching
    crossmatchCatalogues(opt, ref_table, target_table, ref_subtable, target_subtable, ref_finder, target_finder)

if __name__ == '__main__':
    script, ref, target, config = sys.argv
    main(ref, target, config)
