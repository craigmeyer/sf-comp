

class CrossmatchOptions(object):

	def __init__(self, ref, target, filename):

		# Set default values
		self.snap = 0
		self.h = 0.678
		self.boxsize = 0
		self.snapdir = ''
		self.SFdir = ''
		self.SFoutfilename = ''
		self.VELdir = ''
		self.VELoutfilename = ''
		self.catoutfilename = ''
		self.comoving = 1

		with open(filename, 'r') as f:

			for line in f:

				# This line is a comment - skip it
				if line[0] == '#':
					continue

				# Remove all spaces in the line, including leading and trailing whitespace
				line = line.replace(' ', '')
				line = line.strip()

				if not line:
					continue

				# Separate the line into two elements either sign of the assignment operator
				# First element is name of option, second element is the value
				line = line.split('=')

				if line[0] == 'snap':
					self.snap = int(line[1])
				elif line[0] == 'h':
					self.h = float(line[1])
				elif line[0] == 'boxsize':
					self.boxsize = float(line[1])
				elif line[0] == 'snapdir':
					self.snapdir = line[1]
				elif line[0] == 'comoving':
					self.comoving = line[1]
				elif line[0] == 'catoutfilename':
					self.catoutfilename = line[1]

				# SubFind specific options
				elif line[0] == 'SFdir':
					if((ref == 'SF-HBT') or (target == 'SF-HBT')):
						self.SFdir = line[1]
				elif line[0] == 'SFoutfilename':
					if((ref == 'SF-HBT') or (target == 'SF-HBT')):
						self.SFoutfilename = line[1]

				# VELOCIraptor specific options
				elif line[0] == 'VELdir':
					if((ref == 'VEL') or (target == 'VEL')):
						self.VELdir = line[1]
				elif line[0] == 'VELoutfilename':
					if((ref == 'VEL') or (target == 'VEL')):
						self.VELoutfilename = line[1]

				# Handle invalid configuration option
				else:
					raise OSError('Invalid config option passed: %s. Please only use options in the sample config file.' % line[0])





class AnalysisOptions(object):

	def __init__(self, filename):

		# Set default values
		self.snap = 0
		self.h = 0.678
		self.boxsize = 0
		self.snapdir = ''
		self.SFdir = ''
		self.VELdir = ''

		with open(filename, 'r') as f:

			for line in f:

				# This line is a comment - skip it
				if line[0] == '#':
					continue

				# Remove all spaces in the line, including leading and trailing whitespace
				line = line.replace(' ', '')
				line = line.strip()

				if not line:
					continue

				# Separate the line into two elements either sign of the assignment operator
				# First element is name of option, second element is the value
				line = line.split('=')

				if line[0] == 'snap':
					self.snap = int(line[1])
				elif line[0] == 'h':
					self.h = float(line[1])
				elif line[0] == 'boxsize':
					self.boxsize = float(line[1])
				elif line[0] == 'snapdir':
					self.snapdir = line[1]

				# SubFind specific options
				elif line[0] == 'SFdir':
					self.SFdir = line[1]

				# VELOCIraptor specific options
				elif line[0] == 'VELdir':
					self.VELdir = line[1]

				# Handle invalid configuration option
				else:
					raise OSerror('Invalid config option passed: %s. Please only use options in the sample config file.' % line[0])

















