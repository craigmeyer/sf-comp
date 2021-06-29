"""
Analysis of matched halo catalogues from two separate group-finding algorithms

Contains routines to organise and manipulate the data, in addition to
preparing it for diagnostic plots

NOTE: Currently only supports SubFind/SubFind-HBT and VELOCIraptor outputs
NOTE: Currently the substructure hierarchy only goes down to subhaloes

//todo: extend functionality extend to operating with catalogues split across multiple files and multiple snapshots
//todo: extend to HBT group-finding algorithms
//todo: Allow for the subsubhaloes in VELOCIraptor
"""

import sys
import numpy as np
import pandas as pd
import h5py
from ReadConfig import AnalysisOptions
import SubFindIO
import VELOCIraptorIO
import common

# Probably do older routines here (particularly, the calculation of potentials, CoM, etc.)

def getSFpartData(opt, z, halodata, pids, coords, vels, nhalo, nsubhalo):
	"""
	Find particles in each (sub)halo in SubFind halo catalogue

	Paramters:
	opt - Configuration options
	z - Redshift of this snapshot
	halodata - Data from SubFind halo catalogue
	pids, coords, vels - PIDs, coordinates, velocities of all particles in this snapshot
	nhalo, nsubhalo - Number of halos/subhalos in the SubFind halo catalogue

	Returns:
	pids_fof - Array of PIDs for all particles in each (sub)halo
	coords_fof - Array of coordinates for all particles in each (sub)halo
	vels_fof - Array of velocities for all particles in each (sub)halo
	"""

	pids_fof = np.empty(nhalo, dtype = 'object')
	coords_fof = np.empty(nhalo, dtype = 'object')
	vels_fof = np.empty(nhalo, dtype = 'object')

	# Lower and upper indices (first and last particle) for each halo
	linds = halodata['Group/GroupOffsetType'][:,1]
	uinds = linds + halodata['Group/GroupLen']

	# PIDs, coords, vels for all particles in each FOF group (including substructure)
	for ihalo in range(nhalo):
		pids_fof[ihalo] = np.array([part for part in pids[linds[ihalo]:uinds[ihalo]]])
		coords_fof[ihalo] = np.array([coord for coord in coords[linds[ihalo]:uinds[ihalo]]])
		vels_fof[ihalo] = np.array([vel for vel in vels[linds[ihalo]:uinds[ihalo]]])

	# Fix wraparound for haloes located near box edges
	common.fixWraparound(opt, coords_fof, z)

	return pids_fof, coords_fof, vels_fof


def splitSFgroups(halodata, pids_fof, coords_fof, vels_fof, pids, coords, vels, nhalo, nsubhalo):
	"""
	Split SubFind FOF groups into different subsets
	This includes the full FOF group, the main (background/central) subhalo, 
	satellite subhaloes, and diffuse material
	This allows for more flexible comparison between other algorithms,
	which may include different subsets of particles in haloes and/or
	subhaloes
	"""

	firstSub = halodata['Group/GroupFirstSub']
	numSubs = halodata['Group/GroupNsubs']
	group_offset = halodata['Group/GroupOffsetType'][:,1]
	group_npart = halodata['Group/GroupLen']
	sub_npart = halodata['Subhalo/SubhaloLen']
	sub_offset = halodata['Subhalo/SubhaloOffsetType'][:,1]

	central_mask = np.where(firstSub != -1)[0] # FOF groups with a central subhalo
	ncentrals = central_mask.size
	sat_mask = np.where(halodata['Subhalo/SubhaloRankInGr'] > 0)[0]
	nsats = sat_mask.size

	# Arrays to hold different subsets of FOF group particles
	pids_central = np.empty(nhalo, dtype = 'object')
	coords_central = np.empty(nhalo, dtype = 'object')
	vels_central = np.empty(nhalo, dtype = 'object')
	pids_sat = np.empty(nsubhalo, dtype = 'object')
	coords_sat = np.empty(nsubhalo, dtype = 'object')
	vels_sat = np.empty(nsubhalo, dtype = 'object')
	pids_background = np.empty(nhalo, dtype = 'object')
	coords_background = np.empty(nhalo, dtype = 'object')
	vels_background = np.empty(nhalo, dtype = 'object')
	pids_diffuse = np.empty(nhalo, dtype = 'object')
	coords_diffuse = np.empty(nhalo, dtype = 'object')
	vels_diffuse = np.empty(nhalo, dtype = 'object')

	# Split FOF group into central subhalo, diffuse material, and 
	# background (central + diffuse)
	for i in range(ncentrals):
		ihalo = central_mask[i]
		# Central subhalo IDs and coordinates
		id_arr = pids_fof[ihalo][0:sub_npart[firstSub[ihalo]]]
		coord_arr = coords_fof[ihalo][0:sub_npart[firstSub[ihalo]]]
		vel_arr = vels_fof[ihalo][0:sub_npart[firstSub[ihalo]]]
		pids_central[ihalo] = id_arr
		coords_central[ihalo] = coord_arr
		vels_central[ihalo] = vel_arr

		# Now add in diffuse particles (not bound in substructure)
		lastSub = firstSub[ihalo] + numSubs[ihalo] - 1
		id_unbound = np.array([part for part in pids[(sub_offset[lastSub] + sub_npart[lastSub]):
													(group_offset[ihalo] + group_npart[ihalo])]])
		coord_unbound = np.array([coord for coord in coords[(sub_offset[lastSub] + sub_npart[lastSub]):
															(group_offset[ihalo] + group_npart[ihalo])]])
		vel_unbound = np.array([vel for vel in vels[(sub_offset[lastSub] + sub_npart[lastSub]):
													(group_offset[ihalo] + group_npart[ihalo])]])
		# There are diffuse particles in this FOF group
		if id_unbound.size > 0:
			# Populate diffuse particle arrays
			pids_diffuse[ihalo] = id_unbound
			coords_diffuse[ihalo] = coord_unbound
			vels_diffuse[ihalo] = vel_unbound
			# Create background particle (central subhalo + diffuse) arrays
			id_arr = np.concatenate((id_arr, id_unbound))
			coord_arr = np.concatenate((coord_arr, coord_unbound))
			vel_arr = np.concatenate((vel_arr, vel_unbound))
		# Populate background particle arrays
		pids_background[ihalo] = id_arr
		coords_background[ihalo] = coord_arr
		vels_background[ihalo] = vel_arr

	# Split FOF into (satellite) subhaloes
	linds = sub_offset
	uinds = linds + sub_npart
	for i in range(nsats):
		isub = sat_mask[i]
		pids_sat[isub] = np.array([part for part in pids[linds[isub]:uinds[isub]]])
		coords_sat[isub] = np.array([coord for coord in coords[linds[isub]:uinds[isub]]])
		vels_sat[isub] = np.array([vel for vel in vels[linds[isub]:uinds[isub]]])

	# Construct FOF components dictionary
	fofdata = {}
	fofdata['Background/PIDs'] = pids_background
	fofdata['Background/Coordinates'] = coords_background
	fofdata['Background/Velocities'] = vels_background
	fofdata['Central/PIDs'] = pids_central
	fofdata['Central/Coordinates'] = coords_central
	fofdata['Central/Velocities'] = vels_central
	fofdata['Diffuse/PIDs'] = pids_diffuse
	fofdata['Diffuse/Coordinates'] = coords_diffuse
	fofdata['Diffuse/Velocities'] = vels_diffuse
	fofdata['Satellite/PIDs'] = pids_sat
	fofdata['Satellite/Coordinates'] = coords_sat
	fofdata['Satellite/Velocities'] = vels_sat

	return fofdata

def splitVRdata(partdata, halodata, pids_halos, pids, coords, vels, nhalo, nsubhalo):
	"""
	Split VELOCIraptor FOF groups into different subsets
	This includes the full FOF group, the background (field) halo, and
	subhaloes
	This allows for more flexible comparison between other algorithms,
	which may include different subsets of particles in haloes and/or
	subhaloes
	"""

	# Arrays to hold different subsets of FOF group particles
	pids_background = np.array(pids_halos[:nhalo], dtype = 'object')
	coords_background = np.empty(nhalo, dtype = 'object')
	vels_background = np.empty(nhalo, dtype = 'object')
	pids_sub = np.array(pids_halos[nhalo:], dtype = 'object')
	coords_sub = np.empty(nsubhalo, dtype = 'object')
	vels_sub = np.empty(nsubhalo, dtype = 'object')

	# Create analogues to pids_halos (i.e. array where each
	# entry is the coordinates and velocities for each particle
	# in that (sub)halo)
	pids_all = np.concatenate(pids_halos)
	pid_idx = np.argsort(pids)
	pid_sorted = pids[pid_idx]
	match_idx = np.searchsorted(pid_sorted, pids_all)
	idxs = pid_idx[match_idx]
	coords_halos = coords[idxs]
	vels_halos = vels[idxs]

	# Indices that mark the first and last particle in each (sub)halo
	linds = partdata['Offset'] + partdata['Offset_unbound']
	uinds = linds + partdata['Npart']

	# (Field) haloes
	for ihalo in range(nhalo):
		coords_background[ihalo] = np.array([coord for coord in coords_halos[linds[ihalo]:uinds[ihalo]]])
		vels_background[ihalo] = np.array([vel for vel in vels_halos[linds[ihalo]:uinds[ihalo]]])

	# Subhaloes
	for isub in range(nhalo, nhalo + nsubhalo):
		idx = isub - nhalo
		coords_sub[idx] = np.array([coord for coord in coords_halos[linds[isub]:uinds[isub]]])
		vels_sub[idx] = np.array([vel for vel in vels_halos[linds[isub]:uinds[isub]]])

	# Get PIDs of all subhaloes hosted by each field halo to create
	# arrays containing ALL particles in the FOF group
	hostHaloID = halodata['hostHaloID']
	pids_fof = np.empty(nhalo, dtype = 'object')
	coords_fof = np.empty(nhalo, dtype = 'object')
	vels_fof = np.empty(nhalo, dtype = 'object')
	for ihalo in range(nhalo):
		subs = np.where(hostHaloID == ihalo + 1)[0] - nhalo
		if subs.size > 0:
			#print(pids_sub[subs])
			#print(np.concatenate(pids_sub[subs]))
			pids_fof[ihalo] = np.concatenate((pids_background[ihalo], np.concatenate(pids_sub[subs])))
			coords_fof[ihalo] = np.concatenate((coords_background[ihalo], np.concatenate(coords_sub[subs])))
			vels_fof[ihalo] = np.concatenate((vels_background[ihalo], np.concatenate(vels_sub[subs])))
		else: # This halo hosts no subhaloes
			pids_fof[ihalo] = pids_background[ihalo]
			coords_fof[ihalo] = coords_background[ihalo]
			vels_fof[ihalo] = vels_background[ihalo]

	# Construct FOF components dictionary
	fofdata = {}
	fofdata['FOF/PIDs'] = pids_fof
	fofdata['FOF/Coordinates'] = coords_fof
	fofdata['FOF/Velocities'] = vels_fof
	fofdata['Background/PIDs'] = pids_background
	fofdata['Background/Coordinates'] = coords_background
	fofdata['Background/Velocities'] = vels_background
	fofdata['Satellite/PIDs'] = pids_sub
	fofdata['Satellite/Coordinates'] = coords_sub
	fofdata['Satellite/Velocities'] = vels_sub

	return fofdata


def main(config_file):

	# Load in configuration options
	opt = AnalysisOptions(config_file)


	##### SNAPSHOT ######

	# Read in snapshot data
	snap_fn = opt.snapdir + '/snapshot_%03d.hdf5' % opt.snap
	fields = {'PartType1': ('Coordinates', 'ParticleIDs', 'Velocities')}
	snapdata = SubFindIO.readSnapshot(opt, snap_fn, fields)
	pcoords, pids, pvels = snapdata

	# Ensure coordinates and velocities are in standard units
	z = common.getZfromSnap(snap_fn)
	sf_unitinfo = SubFindIO.getUnitInfo(snap_fn)


	##### SubFind #####

	# Read in SubFind halo catalogue (loading all fields)
	halo_fn = opt.SFdir + '/fof_subhalo_tab_%03d.hdf5' % opt.snap
	sf_halodata, sf_nhalo, sf_nsubhalo = SubFindIO.loadSFhaloCatalogue(halo_fn)

	# Get particle data for each SubFind (sub)halo
	sf_pids, sf_coords, sf_vels = getSFpartData(opt, z, sf_halodata, pids, pcoords, pvels, sf_nhalo, sf_nsubhalo)
	sf_coords = common.convLengthUnit(opt, sf_coords, z, sf_unitinfo)
	sf_vels = common.convVelUnit(opt, sf_vels, z, sf_unitinfo)

	# Split each SF FOF group into constituent components
	# These are (from largest to smallest): 
	# FOF group, background (central + diffuse), central, satellites, diffuse
	sf_fof_components = splitSFgroups(sf_halodata, sf_pids, sf_coords, sf_vels, pids, pcoords, pvels, sf_nhalo, sf_nsubhalo)
	sf_fof_components['FOF/PIDs'] = sf_pids
	sf_fof_components['FOF/Coordinates'] = sf_coords
	sf_fof_components['FOF/Velocities'] = sf_vels


	##### VELOCIraptor #####

	# Read in VELOCIraptor halo catalogue (loading all fields)
	base_fn = opt.VELdir + '/snapshot_%03d.VELOCIRAPTOR' % opt.snap
	vr_halodata, _ = VELOCIraptorIO.readVELhalodata(base_fn)
	vr_partdata = VELOCIraptorIO.readVELparticledata(base_fn)
	vr_pids = vr_partdata['Particle_IDs']
	vr_nhalo = np.where(vr_halodata['hostHaloID'] == -1)[0].size
	vr_nsubhalo = vr_halodata['hostHaloID'].size - vr_nhalo

	# Split each VR FOF group into constituent components
	# THese are (from largest to smallest):
	# FOF group, field (central + diffuse) halo, subhaloes
	vr_fof_components = splitVRdata(vr_partdata, vr_halodata, vr_pids, pids, pcoords, pvels, vr_nhalo, vr_nsubhalo)

	# Fix wraparound for VELOCIraptor (sub)halo coordinates
	common.fixWraparound(opt, vr_fof_components['FOF/Coordinates'], z)
	common.fixWraparound(opt, vr_fof_components['Background/Coordinates'], z)
	common.fixWraparound(opt, vr_fof_components['Satellite/Coordinates'], z)

	# Ensure VELOCIraptor coordinates and velocities are in correct units
	vr_unitinfo = VELOCIraptorIO.getUnitInfo(base_fn)
	

	return sf_fof_components, vr_fof_components, sf_halodata, vr_halodata, vr_partdata, sf_unitinfo, vr_unitinfo, sf_pids, vr_pids


if __name__ == '__main__':
	script, config = sys.argv
	main(config)