#!/usr/bin/python

__author__ = 'hofmann'

import unittest
import os
import shutil
import tempfile
from validator import Validator


def touch(file_name, times=None):
	with open(file_name, 'a'):
		os.utime(file_name, times)


class DefaultValidator(unittest.TestCase):
	log_file_path = 'unittest_log.txt'
	taxonomy_directory = '/home/hofmann/Downloads/ncbi-taxonomy_20150130/'
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

if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(TestValidatorMethods)
	unittest.TextTestRunner(verbosity=2, buffer=True).run(suite)
