"""

This script contains plotting routines

//todo: Incorporate a common_plot.py file to have common routines for setting up
plotting figures/windows (e.g. axis, legend, label, tick, etc. setup)
//todo: Document each function
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

def plot_crossmatch_summary(plt, halomatch_cat, submatch_cat, sf_halodata, vr_halodata):
	"""
	Plots ratio of particles in matched (sub)haloes across two catalogues

	NOTE: Only supports SubFind/SubFind-HBT and VELOCIraptor at this time

	//todo: Generalise functionality to other group-finding algorithms
	(plt, halomatch_cat, submatch_cat, g_npart1, g_npart2, s_npart1, s_npart2)
	"""

	fig = plt.figure(figsize=(7,7))

	# //todo: Slice by index rather than column name to generalise
	sf_halos = halomatch_cat['SF-HBT FOF Group']
	vr_halos = halomatch_cat['VEL Field Halo']
	sf_subs = submatch_cat['SF-HBT Subhalo']
	vr_subs = submatch_cat['VEL Subhalo'] + np.where(vr_halodata['hostHaloID'] == -1)[0].size

	group_npart = sf_halodata['Group/GroupLen']
	sub_npart = sf_halodata['Subhalo/SubhaloLen']
	npart = vr_halodata['npart']

	plt.scatter(group_npart[sf_halos], group_npart[sf_halos] / npart[vr_halos], color = 'b', s = 10, label = 'FOF Groups')
	plt.scatter(sub_npart[sf_subs], sub_npart[sf_subs] / npart[vr_subs], color = 'r', s = 10, label = 'Subhaloes')
	plt.hlines(np.mean(group_npart[sf_halos] / npart[vr_halos]), xmin = 0, xmax = 3e5, color = 'b', ls = '--')
	plt.hlines(np.mean(sub_npart[sf_subs] / npart[vr_subs]), xmin = 0, xmax = 3e5, color = 'r', ls = '--')
	plt.xlim(2e1, 2e5)
	plt.ylim(0, 4)
	plt.xlabel(r'$N_{\rm part, SF-HBT}$')
	plt.ylabel(r'$N_{\rm part, SF-HBT}\,/\,N_{\rm part, VEL}$')
	plt.grid(alpha = 0.5)
	plt.xscale('log')
	plt.legend()

	return fig


def plot_viz_halo(plt, coords, index):
	"""
	Visualise a single (sub)halo through three planar projections

	//todo: Add text label for the (sub)halo visualised
	"""

	fig, axes = plt.subplots(1, 3)

	for i in range(axes.size):

		ax = np.ndarray.flatten(axes)[i]
		ax.grid(alpha = 0.5)

		if i == 0:
			ax.scatter(coords[index][:,0], coords[index][:,1], color = 'k', s = 10)
			ax.set_xlabel('X [Mpc / h]')
			ax.set_ylabel('Y [Mpc / h]')
		if i == 1:
			ax.scatter(coords[index][:,0], coords[index][:,2], color = 'k', s = 10)
			ax.set_xlabel('X [Mpc / h]')
			ax.set_ylabel('Z [Mpc / h]')
		if i == 2:
			ax.scatter(coords[index][:,1], coords[index][:,2], color = 'k', s = 10)
			ax.set_xlabel('Y [Mpc / h]')
			ax.set_ylabel('Z [Mpc / h]')

	fig.tight_layout()

	return fig


def plot_offsetPDF(plt, offsets):
	"""
	Plot PDF of offsets between centre estimates between two group-finding
	algorithms

	Examples of metrics include minimum potential, most bound particle,
	and centre-of-mass

	Offsets are in units of a scale radius of the (sub)halo in question 
	(e.g. R200 or half-mass radius)
	"""

	# Construct PDF from kernel density estimate
	kde = gaussian_kde(offsets)
	xarr = np.linspace(-2, 2, 1000)
	plt.plot(xarr, kde(xarr))
	#plt.hist(offsets, bins = 20, range = (-2, 2), density = True, histtype = 'step')

def plot_viz_fof_components(plt, fof_components1, fof_components2, h_ind1, h_ind2, s_inds1, s_inds2, rad1, rad2, cen1, cen2, 
	numsubs = 2):
	"""
	Visualise the different components (whole FOF group, background 
	(main/central) (sub)halo, satellite subhaloes, etc.) for a pair
	of matched FOF groups
	"""

	fig, axes = plt.subplots(3, 2, sharex = True, sharey = True, figsize = (15, 15))


	for i in range(axes.size):

		ax = np.ndarray.flatten(axes)[i]
		ax.grid(alpha = 0.5)
		# Axis labels
		if i > 3: # Bottom row
			xlab = 'X [Mpc / h]'
		else:
			xlab = ''
		ax.set_xlabel(xlab)
		if i % 2 == 0: # Left column
			ax.set_ylabel('Y [Mpc / h]')

		# Background (main/central) halo
		if i > 1:
			if i % 2 == 0:
				ax.scatter(fof_components1['Background/Coordinates'][h_ind1][:,0], 
					fof_components1['Background/Coordinates'][h_ind1][:,1], color = 'b', s = 10, alpha = 0.5)
			else:
				ax.scatter(fof_components2['Background/Coordinates'][h_ind2][:,0], 
					fof_components2['Background/Coordinates'][h_ind2][:,1], color = 'r', s = 10, alpha = 0.5)

		# Plot centre and scale radius (r200c, r200m, r_halfmass, etc.
		if i % 2 == 0:
			ax.scatter(cen1[0], cen1[1], marker = 'x', s = 100, color = 'k', lw = 3)
			circ_rad = plt.Circle((cen1[0], cen1[1]), rad1, edgecolor = 'k', facecolor = 'none', ls = '--', lw = 2)
			ax.add_patch(circ_rad)
		else:
			ax.scatter(cen2[0], cen2[1], marker = 'x', s = 100, color = 'k', lw = 3)
			circ_rad = plt.Circle((cen2[0], cen2[1]), rad2, edgecolor = 'k', facecolor = 'none', ls = '--', lw = 2)
			ax.add_patch(circ_rad)



	axes = np.ndarray.flatten(axes)

	# Entire FOF group (including substructure)
	ax0 = axes[0]
	ax0.scatter(fof_components1['FOF/Coordinates'][h_ind1][:,0], fof_components1['FOF/Coordinates'][h_ind1][:,1],
				color = 'b', s = 10, alpha = 0.5)

	ax1 = axes[1]
	ax1.scatter(fof_components2['FOF/Coordinates'][h_ind2][:,0], fof_components2['FOF/Coordinates'][h_ind2][:,1],
				color = 'r', s = 10, alpha = 0.5)

	# Background (main/central) halo + subhaloes
	ax4 = axes[4]
	if s_inds1.size < numsubs:
		for isub in range(s_inds1.size):
			# Denote the two most massive satellite subhaloes with different colours
			if isub == 0:
				col = 'g'
			elif isub == 1:
				col = 'm'
			else:
				col = 'k'
			ax4.scatter(fof_components1['Satellite/Coordinates'][s_inds1[isub]][:,0], 
				fof_components1['Satellite/Coordinates'][s_inds1[isub]][:,1], color = col, s = 10, alpha = 0.5)
	else:
		for isub in range(numsubs):
			if isub == 0:
				col = 'g'
			elif isub == 1:
				col = 'm'
			else:
				col = 'k'
			ax4.scatter(fof_components1['Satellite/Coordinates'][s_inds1[isub]][:,0], 
				fof_components1['Satellite/Coordinates'][s_inds1[isub]][:,1], color = col, s = 10, alpha = 0.5)

	#//todo: Invesitgate why colouring of most massive subhaloes isn't working with VEL
	ax5 = axes[5]
	if s_inds2.size < numsubs:
		for isub in range(s_inds2.size):
			if isub == 0:
				col = 'g'
			elif isub == 1:
				col = 'm'
			else:
				col = 'k'
			ax5.scatter(fof_components2['Satellite/Coordinates'][s_inds2[isub]][:,0], 
				fof_components2['Satellite/Coordinates'][s_inds2[isub]][:,1], color = col, s = 10, alpha = 0.5)
	else:
		for isub in range(numsubs):
			if isub == 0:
				col = 'g'
			elif isub == 1:
				col = 'm'
			else:
				col = 'k'
			ax5.scatter(fof_components2['Satellite/Coordinates'][s_inds2[isub]][:,0], 
				fof_components2['Satellite/Coordinates'][s_inds2[isub]][:,1], color = col, s = 10, alpha = 0.5)

	# Eliminate extraneous whitespace between panels
	fig.subplots_adjust(hspace = 0, wspace = 0)

	return fig


def plot_fof_pop_props(plt, halomatch_cat, npart, r200c_1, r200c_2, r200m_1, r200m_2, mass1, mass2, merit):

	sf_halos = halomatch_cat['SF-HBT FOF Group']
	vr_halos = halomatch_cat['VEL Field Halo']

	fig, axes = plt.subplots(4, 1, sharex = True, sharey = True, figsize = (15, 15))

	for i in range(axes.size):

		# Configure axes
		ax = np.ndarray.flatten(axes)[i]
		ax.set_xscale('log')
		ax.set_xlim(2e1, 2e5)
		ax.set_ylim(0, 2)
		ax.grid(axis = 'y', alpha = 0.5)
		if i == 3:
			xlab = r'$N_{\rm part, SF}$'
		else:
			xlab = ''
		ax.set_xlabel(xlab)
		# 100 particle marker
		ax.axvline(x = 100, ymin = 0, ymax = 2, color = 'g', lw = 2, ls = '--')

	axes = np.ndarray.flatten(axes)

	#//todo: Test if fontsize needs to be adjusted
	ax0 = axes[0]
	ax0.scatter(npart[sf_halos], r200c_1[sf_halos] / r200c_2[vr_halos], color = 'k')
	ax0.set_ylabel(r'$R_{\rm 200c, SF} / R_{\rm 200c, VEL}$')

	ax1 = axes[1]
	ax1.scatter(npart[sf_halos], r200m_1[sf_halos] / r200m_2[vr_halos], color = 'k')
	ax1.set_ylabel(r'$R_{\rm 200m, SF} / R_{\rm 200m, VEL}$')

	ax2 = axes[2]
	ax2.scatter(npart[sf_halos], mass1[sf_halos] / mass2[sf_halos], color = 'k')
	ax2.set_ylabel(r'$M_{\rm SF} / M_{\rm VEL}$')

	ax3 = axes[3]
	ax3.scatter(npart[sf_halos], merit, color = 'k')
	ax3.set_ylabel('Merit Score')

	fig.subplots_adjust(hspace = 0, wspace = 0)

	return fig


def plot_subhalo_pop_props(plt, submatch_cat, npart, vmax1, vmax2, rhalf1, rhalf2, mass1, mass2, merit):

	sf_subs = submatch_cat['SF-HBT Subhalo']
	vr_subs = submatch_cat['VEL Subhalo'] + np.where(vr_halodata['hostHaloID'] == -1)[0].size

	fig, axes = plt.subplots(4, 1, sharex = True, sharey = True)

	for i in range(axes.size):

		# Configure axes
		ax = np.ndarray.flatten(axes)[i]
		ax.set_xscale('log')
		ax.set_xlim(2e1, 2e5)
		ax.set_ylim(0, 2)
		ax.grid(axis = 'y', alpha = 0.5)
		if i == 3:
			xlab = r'$N_{\rm part, SF}$'
		else:
			xlab = ''
		ax.set_xlabel(xlab)
		# 100 particle marker
		ax.axvline(x = 100, ymin = 0, ymax = 2, color = 'g', lw = 2, ls = '--')

	axes = np.ndarray.flatten(axes)

	#//todo: Test if fontsize needs to be adjusted
	ax0 = axes[0]
	ax0.scatter(npart[sf_subs], vmax1[sf_subs] / vmax2[vr_subs], color = 'k')
	ax0.set_ylabel(r'$V_{\rm max, SF} / V_{\rm max, VEL}$')

	ax1 = axes[1]
	ax1.scatter(npart[sf_subs], rhalf1[sf_subs] / rhalf2[vr_subs], color = 'k')
	ax1.set_ylabel(r'$R_{\rm 1/2, SF} / R_{\rm 1/2, VEL}$')

	ax2 = axes[2]
	ax2.scatter(npart[sf_subs], mass1[sf_subs] / mass2[sf_subs], color = 'k')
	ax2.set_ylabel(r'$M_{\rm SF} / M_{\rm VEL}$')

	ax3 = axes[3]
	ax3.scatter(npart[sf_subs], merit, color = 'k')
	ax3.set_ylabel('Merit Score')

	fig.subplots_adjust(hspace = 0, wspace = 0)


def plot_inout_PDFs(plt, npart_in1, npart_in2, npart_out1, npart_out2, merit_in, merit_out):
	"""
	Plot PDFs comparing the particle number and merit score
	in and inner and outer region

	//todo: Add legends
	//todo: Maybe configure so that inner and outer radii are parameters
	and construct arrays within this function
	//todo: See if histograms look better than KDE PDF
	"""

	fig, axes = plt.subplots(1, 2, sharey = True)

	for i in range(axes.size):

		# Configure axes
		ax = np.ndarray.flatten(axes)[i]
		ax.grid(alpha = 0.5)
		if i == 0:
			ylab = 'Density'
		else:
			ylab = ''
		ax.set_ylabel(ylab)

		if i == 0:
			# Construct PDFs from kernel density estimate
			kde_in = gaussian_kde(npart_in1 / npart_in2)
			kde_out = gaussian_kde(npart_out1 / npart_out2)
			xarr = np.linspace(0, 2, 1000)
			ax.plot(xarr, kde_in(xarr), color = 'b')
			ax.plot(xarr, kde_out(xarr), color = 'r')
			#ax.hist(npart_in1 / npart_in2, bins = 20, range = (0, 2), density = True, histtype = 'step')
			#ax.hist(npart_out1 / npart_out2, bins = 20, range = (0, 2), density = True, histtype = 'step')
		else:
			kde_in = gaussian_kde(merit_in)
			kde_out = gaussian_kde(merit_out)
			xarr = np.linspace(0, 1, 1000)
			ax.plot(xarr, kde_in(xarr), color = 'b')
			ax.plot(xarr, kde_out(xarr), color = 'r')
			#ax.hist(merit_in, bins = 10, range = (0, 1), density = True, histtype = 'step')
			#ax.hist(merit_out, bins = 10, range = (0, 1), density = True, histtype = 'step')

	fig.subplots_adjust(hspace = 0, wspace = 0)


def plot_hmf(plt):

	print('dummy function')

def plot_cmf(plt):

	print('dummy function')




