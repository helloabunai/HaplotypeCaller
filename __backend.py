import os

class Colour:

	def __init__(self):
		pass

	purple = '\033[95m'
	cyan = '\033[96m'
	darkcyan = '\033[36m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	bold = '\033[1m'
	underline = '\033[4m'
	end = '\033[0m'

class HaploCalcException:
	def __init__(self):
		pass

class IndividualSNP:
	def __init__(self):
		self.SNPID = ''
		self.present_haplotypes = []
		self.haplotype_frequency = []

	## setters
	def set_SNPID(self, value):
		self.SNPID = value
	def append_haplotype(self, value):
		self.present_haplotypes.append(value)
	def set_frequency(self, value):
		self.haplotype_frequency = value


	## getters
	def get_SNPID(self):
		return self.SNPID
	def get_haplotypes(self):
		return self.present_haplotypes
	def get_frequency(self):
		return self.haplotype_frequency

class IndividualSubject:
	def __init__(self):
		self.subjectID = ''
		self.age = ''
		self.sex = ''
		self.mutations = {}

	## setters
	def set_subjectID(self, value):
		self.subjectID = value
	def set_age(self, value):
		self.age = value
	def set_sex(self, value):
		self.sex = value
	def append_mutation(self, snp_key, mutation_value):
		try:
			self.mutations[snp_key].append(mutation_value)
		except KeyError:
			self.mutations[snp_key] = [mutation_value]

	## getters
	def get_subjectID(self):
		return self.subjectID
	def get_subjectAge(self):
		return self.age
	def get_subjectSex(self):
		return self.sex
	def get_mutations(self):
		return self.mutations


def sanitise_io(candidate_input, candidate_output):

	## input file
	if not os.path.isfile(candidate_input):
		raise HaploCalcException("Specified input path does not exist!")
	if not candidate_input.endswith(".csv"):
		raise HaploCalcException("Specified input file is not a CSV!")
	if not os.path.exists(candidate_output):
		force_mkdir(candidate_output)
		return True

	return True

def force_mkdir(path):
	try:
		os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise
