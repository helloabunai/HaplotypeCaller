#/usr/bin/python
__version__ = 0.1
__author__ = 'alastair.maxwell@glasgow.ac.uk'

##
## Python libraries
import os
import sys
import pandas
import argparse
import collections
import logging as log
import multiprocessing

## global threads
THREADS = multiprocessing.cpu_count()

##
## Backend junk
from .__backend import sanitise_io
from .__backend import Colour as clr
from .__backend import IndividualSNP
from .__backend import IndividualSubject
from .__backend import HaploCalcException

class HaplotypeCalculator:
	def __init__(self):
		
		"""
		Docstring __todo__ lmAO
		"""

		##
		## Argument parser from CLI
		self.parser = argparse.ArgumentParser(prog='haplocalc', description='HaploCalc. Calculates most common haplotypes from SNP data.')
		self.parser.add_argument('-v', '--verbose', help='Verbose output mode. Setting this flag enables verbose output. Default: off.', action='store_true')
		self.parser.add_argument('-i', '--input', help='Input CSV file of SNP data for calculation', nargs=1, required=True)
		self.parser.add_argument('-t', '--threads', help='Thread utilisation. Each sample present in the input will be sent to a different processor, if specified.', type=int, choices=xrange(1, THREADS+1), default=THREADS)
		self.parser.add_argument('-o', '--output', help='Output path. Specify a directory you wish output to be directed towards.', metavar='output', nargs=1, required=True)
		self.args = self.parser.parse_args()
		self.header = ''

		##
		## Set verbosity for CLI output
		if self.args.verbose:
			log.basicConfig(format='%(message)s', level=log.DEBUG)
			log.info('{}{}{}{}'.format(clr.bold, 'hplc__ ', clr.end, 'HaploCalc: Calculates most common haplotypes within SNP data.'))
			log.info('{}{}{}{}'.format(clr.bold, 'hplc__ ', clr.end, 'alastair.maxwell@glasgow.ac.uk\n'))
		else:
			log.basicConfig(format='%(message)s')

		## instance i/o
		self.input_file = self.args.input[0]
		self.output_path = self.args.output[0]

		## check input_file exists, make outdir if doesn't exist
		if sanitise_io(self.input_file, self.output_path):
			parsed_inputname = self.input_file.split('/')[-1].split('.')[0]
			self.mutation_output = os.path.join(self.output_path, '{}_HaploCalcResults.csv'.format(parsed_inputname))
			self.subject_output = os.path.join(self.output_path, '{}_HaploCalcResults_PerSubject.csv'.format(parsed_inputname))

		## dict to store SNP objects, which will contain scraped data
		## also dict for subject objects
		self.instance_mutation_population = {}
		self.instance_subject_population = {}

		## workflow for HaploCalc begins here
		## check our input file has correct column headers
		self.mutation_dataframe = self.validate_headers(self.input_file)
		self.clean_data()
		## organise data by SNP-basis
		self.assign_SNP()
		self.count_haplotypes()
		## organise data by subject-basis
		self.assign_subject()
		self.process_output()
		log.info('{}{}{}{}'.format(clr.green,'hplc__ ',clr.end,'Completed workflow! Exiting..'))

	def validate_headers(self, infi):

		"""
		Check our input file for the required headers in each column for processing
		if valid, then we return the pandas dataframe so it can be processed by other methods
		returns: pd.dataframe(infi)
		"""

		## inform
		log.info('{}{}{}{}'.format(clr.yellow, 'hplc__ ', clr.end, 'Reading input...'))

		## Attempt creation of dataframe
		candidate_dataframe = pandas.read_csv(infi)
		## scrape headers present in each column, create list of expected values
		column_headers = candidate_dataframe.columns.values.tolist()
		required_headers = ['DaughterPlate', 'MasterPlate', 'MasterWell', 'Call', 'X', 'Y', 'SNPID',
		 'SubjectID', 'Age', 'Sex', 'CTG repeat allele1 allele2', 'Unnamed: 11']

		## iterate over both lists at once, comparing each value in step
		## if they match, continue to next value
		## if they don't match, raise an exception saying why
		for value, standard in zip(column_headers, required_headers):
			if value == standard: pass
			else: raise HaploCalcException('Column header mismatch found: Got {}, expected {}.'.format(value, standard))

		## inform and return, if we reached this far
		log.info('{}{}{}{}'.format(clr.green, 'hplc__ ', clr.end, 'Parsed input data successfully!'))
		return candidate_dataframe
		
	def clean_data(self):

		## loop over every row in our dataframe object
		## remove entries where SNP data is missing for whatever well in the plate
		## that this entry represents
		for index, row in self.mutation_dataframe.iterrows():
			if pandas.isnull(row['SubjectID']):
				self.mutation_dataframe.drop(index, inplace=True)

	def assign_SNP(self):

		"""
		Method to iterate over our dataframe, which now contains only full-entries
		Assign values for the current SNP we're working on into objects for later processing
		"""

		log.info('{}{}{}{}'.format(clr.green,'hplc__ ',clr.end,'Assigning data on SNP-basis..'))

		for index, row in self.mutation_dataframe.iterrows():
			current_SNP = row['SNPID']
			current_haplotype = row['CTG repeat allele1 allele2']

			## some entries have slashes instead of spaces, for some fucking reason
			if '/' in current_haplotype:
				current_haplotype = current_haplotype.replace('/', ' ')

			## if the current SNP already has an object made in our instance dictionary of SNPs
			## append the current haplotype for this new observation
			## otherwise, no object has been made; make one and append observation data
			if current_SNP in self.instance_mutation_population:
				self.instance_mutation_population[current_SNP].append_haplotype(current_haplotype)
			else:
				snp_object = IndividualSNP()
				snp_object.set_SNPID(current_SNP)
				snp_object.append_haplotype(current_haplotype)
				self.instance_mutation_population[current_SNP] = snp_object

	def count_haplotypes(self):

		"""
		self explanatory. count haplotypes.
		"""

		for current_SNP, mutation_object in self.instance_mutation_population.iteritems():
			present_haplotypes = mutation_object.get_haplotypes()
			frequency = collections.Counter(present_haplotypes).most_common()
			mutation_object.set_frequency(frequency)

	def assign_subject(self):

		"""
		Method to iterate over our dataframe, which now contains only full-entries
		Assign values for the current SUBJECT we're working on into objects for later processing
		"""

		log.info('{}{}{}{}'.format(clr.green,'hplc__ ',clr.end,'Assigning data on Subject-basis..'))

		for index, row in self.mutation_dataframe.iterrows():
			current_subject = row['SubjectID']
			current_age = row['Age']
			current_sex = row['Sex']
			current_SNP = row['SNPID']
			current_haplotype = row['CTG repeat allele1 allele2']

			## some entries have slashes instead of spaces, for some fucking reason
			if '/' in current_haplotype:
				current_haplotype = current_haplotype.replace('/', ' ')

			if current_subject in self.instance_subject_population:
				self.instance_subject_population[current_subject].append_mutation(current_SNP, current_haplotype)
			else:
				subject_object = IndividualSubject()
				subject_object.set_subjectID(current_subject)
				subject_object.set_age(current_age)
				subject_object.set_sex(current_sex)
				subject_object.append_mutation(current_SNP, current_haplotype)
				self.instance_subject_population[current_subject] = subject_object

	def process_output(self):
		
		## inform
		log.info('{}{}{}{}'.format(clr.green, 'hplc__ ', clr.end, 'Writing output..'))

		## Output our SNP-basis data
		with open(self.mutation_output, 'w') as mutation_outfi:
			mutation_outfi.write('{}, {}, {}\n'.format('SNP ID', 'Haplotype', 'Frequency'))
			for snp_key, mutation_value in self.instance_mutation_population.iteritems():
				snp = mutation_value.get_SNPID()
				frequency = mutation_value.get_frequency()
				for haplotype_frequency in frequency:
					mutation_outfi.write('{}, {}, {}\n'.format(snp, haplotype_frequency[0], haplotype_frequency[1]))
		## Output our Subject-basis data
		with open(self.subject_output, 'w') as subject_outfi:
			subject_outfi.write('{}, {}, {}, {}, {}\n'.format('Subject ID', 'Age', 'Sex', 'SNP ID', 'Haplotype'))
			for subject_key, subject_value in self.instance_subject_population.iteritems():
				subject = subject_value.get_subjectID()
				age = subject_value.get_subjectAge()
				sex = subject_value.get_subjectSex()
				mutations = subject_value.get_mutations()
				for snp_k, mut_v in mutations.iteritems():
					mut_str = ''.join(mut_v)
					subject_outfi.write('{}, {}, {}, {}, {}\n'.format(subject, age, sex, snp_k, mut_str))


def main():
	try:
		HaplotypeCalculator()
	except KeyboardInterrupt:
		log.error('{}{}{}{}'.format(clr.red,'hplc__ ',clr.end,'Fatal: Keyboard Interrupt detected. Exiting.'))
		sys.exit(2)