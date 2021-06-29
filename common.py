"""
Common routines used regularly for analysis

"""

import numpy as np
import h5py

# Constants for unit conversion
CM_TO_MPC = 3.085678e24
CMPERS_TO_KMPERS = 1e5
MSUN_TO_GRAM = 1.989e33

def convLengthUnit(opt, data, z, unitinfo, c_bool = True, h_bool = True):

	length_unit = unitinfo['UnitLength_in_cm']
	data *= length_unit / CM_TO_MPC # Convert to Mpc

	# Units don't include little h scaling
	# Convert from [Mpc] to [Mpc / h]
	if not h_bool:
		data *= opt.h

	# Convert to comoving units, if required
	if not c_bool:
		data /= 1 / (1 + z)

	return data


def convMassUnit(opt, data, unitinfo, h_bool = True):

	mass_unit = unitinfo['UnitMass_in_g']
	data *= mass_unit / MSUN_TO_GRAM / 1e10 # Convert to 10^10 Msun

	# Units don't include little h scaling
	# Convert from [10^10 Msun] to [10^10 Msun/h]
	if not h_bool:
		data *= opt.h

	return data


def convVelUnit(opt, data, z, unitinfo, c_bool = True):

	vel_unit = unitinfo['UnitVelocity_in_cm_per_s']
	data *= vel_unit / CMPERS_TO_KMPERS # Convert to km/s

	# Convert to comoving units, if required
	if not c_bool:
		data /= (1 / (1 + z))**2

	return data


def fixWraparound(opt, coord_arr, z):
	"""
	Corrects particle coordinates for artifical wraparound due to
	being located near the edges of the simulation box

	opt - Configuration options
	coord_arr - Array of coordinates - each entry contains the 
	coordinates of all particles in a given (sub)halo
	"""

	# Scale factor - needed to account for changing box-size at snapshots other
	# than z = 0
	a = 1 / (1 + z)

	# Define wraparound if the x-, y-, or z-, coordinates in a given (sub)halo
	# have a range at least 90% the box size
	threshold = 0.9 * opt.boxsize

	for coords in coord_arr:
		diffs = np.ndarray.max(coords, axis = 0) - np.ndarray.min(coords, axis = 0)
		max_diff = np.max(diffs)

		if max_diff > threshold:
			max_ind = np.where(diffs == max_diff)[0]
			# Determine whether wraparound occurred at low- or high- edge of box
			# and coorect accordingly
			n_largevals = np.where(coords[:, max_ind] > threshold)[0].size
			n_smallvals = np.where(coords[:, max_ind] < threshold)[0].size
			if n_largevals > n_smallvals:
				coords[:, max_ind] = np.maximum(np.abs(coords[:, max_ind] - opt.boxsize * a), coords[:, max_ind])
			else:
				coords[:, max_ind] = np.minimum(np.abs(coords[:, max_ind] - opt.boxsize * a), coords[:, max_ind])


def calcMeritScore(ref_table, target_table, num_matches, ref_bg_pids, target_bg_pids, ref_fof_pids):
	"""
	Calculates the merit score between two matched (sub)haloes according to
	the metric given by Eq. 5 from Poulton et al. 2019

	Parameters:
	ref_table - Reference halo catalogue table of matched objects
	target_table - Target halo catalogue table of mathced objects
	num_matches - Number of matches between the two halo catalogues
	ref_pids - PIDs for each (sub)halo in reference catalogue
	target_pids - PIDs for each (sub)halo in target catalogue

	Returns:
	merit - Array containing the merit scores for each match

	//todo: Maybe extend so can calculate merit score within certain regions
	of the (sub)halo (i.e. r < Rin or r > Rout)
	"""

	merit = np.zeros(num_matches, dtype = 'float64')

	for i in range(num_matches):
		ref = ref_table[i]
		target = target_table[i]
		if ref_bg_pids[ref] is None: # SF halo has no central subhalo, FOF group is entire hierarchy
			p_ref = ref_fof_pids[ref]
		else:
			p_ref = ref_bg_pids[ref]
		p_target = target_bg_pids[target]
		n_ref = p_ref.size
		n_target = p_target.size
		n_common = np.sum(np.in1d(p_ref, p_target))
		merit[i] = n_common**2 / (n_ref * n_target)

	return merit


def getInteriorInds(cen, coords, rin):
	"""
	Get indices of particles in (sub)halo located within r < rin

	Parameters:
	cen - Centre (minimum potential) of (sub)halo
	coords - Coordinates of all particles in the (sub)halo

	Returns:
	in_inds - Indices of all particles which satisfy r < rin
	"""

	dr = cen - coords
	dists = np.sqrt(np.sum(dr**2, axis = 1))
	in_inds = np.where(dists < rin)[0]

	return in_inds


def getExteriorInds(cen, coords, rout):
	"""
	Get indices of particles in (sub)halo located within r > rout

	Parameters:
	cen - Centre (minimum potential) of (sub)halo
	coords - Coordinates of all particles in the (sub)halo

	Returns:
	out_inds - Indices of all particles which satisfy r > rout
	"""

	dr = cen - coords
	dists = np.sqrt(np.sum(dr**2, axis = 1))
	out_inds = np.where(dists > rin)[0]

	return out_inds


def getZfromSnap(fn):

	snapfile = h5py.File(fn, 'r')

	header = snapfile['Header'].attrs
	z = header['Redshift']
	snapfile.close()

	return z

def calc_coord_offsets(coords1, coords2):
	""" 
	Determine distance between two sets of coordaintes.
	Coordinates passed can be two separate coordinates, or two sets of many coordinates, 
	in which case the distance is computed between the corresponding entries in the two arrays

	Parameters:
	coords1 - First coordinate (or array of coordinates)
	coords2 - Second coordinate (or array of coordinates)

	Returns:
	Distance between the two coordinates, or array of distances between many coordinates
	"""

	if len(coords1.shape) == 1:
		return np.linalg.norm(coords1 - coords2)
	else:
		return np.linalg.norm(coords1 - coords2, axis = 1)
















