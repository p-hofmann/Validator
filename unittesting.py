#!/usr/bin/python

__author__ = 'hofmann'

import unittest
import os
import sys
import shutil
import tempfile
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from validator import Validator
from sequencevalidator import SequenceValidator


def touch(file_name, times=None):
	with open(file_name, 'a'):
		os.utime(file_name, times)


class DefaultValidator(unittest.TestCase):
	log_file_path = 'unittest_log.txt'
	file_names = ["f0", "f1", "f2"]
	executables = ["unittesting.py"]
	valid_numbers = [4, 0, -1, 1.2, 999.222, 9999999999999999999.09999]
	invalid_numbers = ["4", "", None]

	def setUp(self):
		self.file_stream = open(DefaultValidator.log_file_path, 'a')
		self.validator = Validator(logfile=self.file_stream, verbose=False)
		self.temp_directory = tempfile.gettempdir()
		self.directory = tempfile.mkdtemp(dir=self.temp_directory)
		self.sub_directory = tempfile.mkdtemp(dir=self.directory)
		for file_name in DefaultValidator.file_names:
			file_path = os.path.join(self.directory, file_name)
			touch(file_path)

	def tearDown(self):
		self.validator.close()
		self.validator = None
		self.file_stream.close()
		self.file_stream = None
		if os.path.isfile(DefaultValidator.log_file_path):
			os.remove(DefaultValidator.log_file_path)
		if os.path.isdir(self.directory):
			shutil.rmtree(self.directory)


class TestValidatorMethods(DefaultValidator):

	def test_file_validation(self):
		# test for existence
		for file_name in DefaultValidator.file_names:
			file_path = os.path.join(self.directory, file_name)
			self.assertTrue(
				self.validator.validate_file(file_path, executable=False, silent=False),
				"File existence was not confirmed")
			self.assertFalse(
				self.validator.validate_file(file_path, executable=True, silent=True),
				"Non executable file was wrongly confirmed as executable")
		self.assertFalse(
			self.validator.validate_file("made_up_name", executable=False, silent=True),
			"Non existence file was confirmed")

		# test for executable right
		for file_name in DefaultValidator.executables:
			self.assertTrue(
				self.validator.validate_file(file_name, executable=True, silent=False),
				"Executable was not confirmed as such")

	def test_directory_validation(self):
			self.assertTrue(
				self.validator.validate_dir(
					self.directory, only_parent=True, sub_directories=None, file_names=None, key=None, silent=False),
				"ParentDir existence was not confirmed")
			self.assertTrue(
				self.validator.validate_dir(
					self.directory, only_parent=False, sub_directories=None, file_names=None, key=None, silent=False),
				"Dir existence was not confirmed")
			sub_directories = [os.path.basename(self.sub_directory)]
			self.assertTrue(
				self.validator.validate_dir(
					self.directory, only_parent=False, sub_directories=sub_directories),
				"SubDir existence was not confirmed")
			self.assertTrue(
				self.validator.validate_dir(
					self.directory, only_parent=False, file_names=self.file_names),
				"Existence of files in dir was not confirmed")

	def test_full_path(self):
		paths = [
			(".", os.getcwd()),
			("./", os.getcwd()),
			("~", "/home/hofmann"),
			("~/", "/home/hofmann"),
			("/afg/", "/afg"),
			("/afg/qwe", "/afg/qwe"),
			("qwe", os.path.join(os.getcwd(), "qwe")),
		]
		for path in paths:
			full_path = self.validator.get_full_path(path[0])
			self.assertEqual(full_path, path[1], "'{}' != '{}'".format(full_path, path[1]))

	def test_number_validation(self):
		for value in TestValidatorMethods.valid_numbers:
			self.assertTrue(
				self.validator.validate_number(value, minimum=None, maximum=None, zero=True, key=None, silent=False)
				)

		for value in TestValidatorMethods.invalid_numbers:
			with self.assertRaises(AssertionError):
				self.validator.validate_number(value, minimum=None, maximum=None, zero=True, key=None, silent=False)

		self.assertTrue(
			self.validator.validate_number(0, minimum=None, maximum=None, zero=True, key=None, silent=False)
			)
		self.assertFalse(
			self.validator.validate_number(0, minimum=None, maximum=None, zero=False, key=None, silent=True)
			)
		self.assertFalse(
			self.validator.validate_number(0, minimum=1, maximum=None, zero=False, key=None, silent=True)
			)
		self.assertFalse(
			self.validator.validate_number(0, minimum=None, maximum=-1, zero=False, key=None, silent=True)
			)

	def test_free_space_validation(self):
		self.assertTrue(
			self.validator.validate_free_space(
				self.directory,
				required_space_in_bytes=0,
				silent=False)
			)

		self.assertTrue(
			self.validator.validate_free_space(
				self.directory,
				required_space_in_kb=0,
				silent=False)
			)

		self.assertTrue(
			self.validator.validate_free_space(
				self.directory,
				required_space_in_mb=0,
				silent=False)
			)

		self.assertTrue(
			self.validator.validate_free_space(
				self.directory,
				required_space_in_gb=0,
				silent=False)
			)

		with self.assertRaises(AssertionError):
			self.assertTrue(
				self.validator.validate_free_space(
					self.directory,
					required_space_in_bytes=0, required_space_in_kb=0,
					silent=False)
				)

		free_space = self.validator.free_space_in_giga_bytes(self.directory)
		self.assertFalse(
			self.validator.validate_free_space(
				self.directory,
				required_space_in_gb=free_space+1,
				silent=True),
			"Wrongly {} > {} as True declared".format(free_space, free_space+1)
			)


class DefaultSequenceValidator(unittest.TestCase):
	_test_case_id = 0
	_success = False

	dir_input = "unittest_input"
	dir_output = "unittest_output_fa_{}"

	log_filename = 'unittest_log.txt'
	filename_fasta = "unittest.fasta"
	filename_fasta_bad0 = "bad0.fasta"
	filename_fasta_bad1 = "bad1.fasta"
	filename_fasta_ambiguous = "unittest_ambiguous.fasta"
	filename_fq = "unittest.fq"

	def setUp(self):
		self.dir_output = self.dir_output.format(DefaultSequenceValidator._test_case_id)
		if os.path.isdir(self.dir_output):
			shutil.rmtree(self.dir_output)
		os.mkdir(self.dir_output)
		sys.stderr.write("\n{}... ".format(DefaultSequenceValidator._test_case_id)),
		DefaultSequenceValidator._test_case_id += 1

		logfile = os.path.join(self.dir_output, self.log_filename)
		self.file_stream = open(logfile, 'a')
		self.validator = SequenceValidator(logfile=self.file_stream, verbose=False)

	def tearDown(self):
		self.validator.close()
		self.validator = None
		self.file_stream.close()
		self.file_stream = None
		if self._success:
			shutil.rmtree(self.dir_output)


class TestSequenceValidatorMethods(DefaultSequenceValidator):

	def test_file_validation(self):
		fasta_file = os.path.join(self.dir_input, self.filename_fasta)
		fasta_file_ambiguous = os.path.join(self.dir_input, self.filename_fasta_ambiguous)
		fasta_file_bad0 = os.path.join(self.dir_input, self.filename_fasta_bad0)
		fasta_file_bad1 = os.path.join(self.dir_input, self.filename_fasta_bad1)
		fastq_file = os.path.join(self.dir_input, self.filename_fq)
		self.assertTrue(self.validator.validate_sequence_file(fasta_file, "fasta", "dna", ambiguous=False))
		self.assertFalse(self.validator.validate_sequence_file(fasta_file_ambiguous, "fasta", "dna", ambiguous=False, silent=True))
		self.assertTrue(self.validator.validate_sequence_file(fasta_file_ambiguous, "fasta", "dna", ambiguous=True))
		self.assertFalse(self.validator.validate_sequence_file(fasta_file_bad0, "fasta", "dna", ambiguous=False, silent=True))
		self.assertFalse(self.validator.validate_sequence_file(fasta_file_bad1, "fasta", "dna", ambiguous=False, silent=True))
		self.assertTrue(self.validator.validate_sequence_file(fastq_file, "fastq", "dna", ambiguous=False))
		self._success = True

	def test_sequence_validation(self):
		fasta_file = os.path.join(self.dir_input, self.filename_fasta)
		fasta_file_ambiguous = os.path.join(self.dir_input, self.filename_fasta_ambiguous)
		fastq_file = os.path.join(self.dir_input, self.filename_fq)

		# IUPAC.ambiguous_rna
		with open(fasta_file) as stream:
			for seq_record in SeqIO.parse(stream, "fasta", alphabet=IUPAC.unambiguous_dna):
				self.assertTrue(self.validator.validate_sequence(seq_record.seq))

		with open(fasta_file_ambiguous) as stream:
			for seq_record in SeqIO.parse(stream, "fasta", alphabet=IUPAC.ambiguous_dna):
				self.assertTrue(self.validator.validate_sequence(seq_record.seq))

		with open(fastq_file) as stream:
			for seq_record in SeqIO.parse(stream, "fastq", alphabet=IUPAC.unambiguous_dna):
				self.assertTrue(self.validator.validate_sequence(seq_record.seq))
		self._success = True

	def test_sequence_id_validation(self):
		fasta_file = os.path.join(self.dir_input, self.filename_fasta)
		fastq_file = os.path.join(self.dir_input, self.filename_fq)

		# IUPAC.ambiguous_rna
		with open(fasta_file) as stream:
			for seq_record in SeqIO.parse(stream, "fasta", alphabet=IUPAC.unambiguous_dna):
				self.assertTrue(self.validator.validate_sequence_id(seq_record.id))

		with open(fastq_file) as stream:
			for seq_record in SeqIO.parse(stream, "fastq", alphabet=IUPAC.unambiguous_dna):
				self.assertTrue(self.validator.validate_sequence_id(seq_record.id))
		self._success = True

	def test_sequence_description_validation(self):
		fasta_file = os.path.join(self.dir_input, self.filename_fasta)
		fastq_file = os.path.join(self.dir_input, self.filename_fq)

		# IUPAC.ambiguous_rna
		with open(fasta_file) as stream:
			for seq_record in SeqIO.parse(stream, "fasta", alphabet=IUPAC.unambiguous_dna):
				self.assertTrue(self.validator.validate_sequence_description(seq_record.description))

		with open(fastq_file) as stream:
			for seq_record in SeqIO.parse(stream, "fastq", alphabet=IUPAC.unambiguous_dna):
				self.assertTrue(self.validator.validate_sequence_description(seq_record.description))
		self._success = True

if __name__ == '__main__':
	suite0 = unittest.TestLoader().loadTestsFromTestCase(TestValidatorMethods)
	suite1 = unittest.TestLoader().loadTestsFromTestCase(TestSequenceValidatorMethods)
	alltests = unittest.TestSuite([suite0, suite1])
	unittest.TextTestRunner(verbosity=2).run(alltests)
