__author__ = 'hofmann'
__version__ = '0.0.1'

import os
import math
from numbers import Number
from scripts.loggingwrapper import LoggingWrapper


class Validator(object):

	_map_logfile_handler = dict()

	def __init__(self, logfile=None, verbose=False):
		"""
			Collection of methods for value validations

			@attention: config_file argument may be file path or stream.

			@param logfile: file handler or file path to a log file
			@type logfile: file | FileIO | None
			@param verbose: Not verbose means that only warnings and errors will be past to stream
			@type verbose: bool

			@return: None
			@rtype: None
		"""
		self._logger = LoggingWrapper("Validator", verbose=verbose)
		if logfile:
			self._logger.set_log_file(logfile)

	def __exit__(self, type, value, traceback):
		self.close()

	def __enter__(self):
		return self

	def close(self):
		self._logger.close()
		self._logger = None

	def validate_file(self, file_path, executable=False, key=None, silent=False):
		"""
			Collection of methods for value validations

			@attention: config_file argument may be file path or stream.

			@param file_path: path to a file
			@type file_path: basestring
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		assert isinstance(executable, bool)
		assert isinstance(silent, bool)
		assert key is None or isinstance(key, basestring)
		assert file_path is None or isinstance(file_path, basestring)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if file_path is None:
			if not silent:
				self._logger.error("{}Invalid directory".format(prefix))
			return False

		file_path = self.get_full_path(file_path)
		parent_directory = os.path.dirname(file_path)
		if not self.validate_dir(parent_directory, key=key):
			return False

		if not os.path.isfile(file_path):
			if not silent:
				self._logger.error("{}File does not exist: '{}'".format(prefix, file_path))
			return False

		if executable and not os.access(file_path, os.X_OK):
			if not silent:
				self._logger.error("{}Permission error. File can not be executed '{}'".format(prefix, file_path))
			return False
		return True

	def validate_dir(self, directory, only_parent=False, sub_directories=None, file_names=None, key=None, silent=False):
		"""
			Validate existence of directory or parent directory or sub directories and files.

			@attention:

			@param directory: directory path of a folder
			@type directory: basestring
			@param only_parent: test only the existence of the parent directory
			@type only_parent: bool
			@param sub_directories: test the existence of sub directories
			@type sub_directories: list[basestring]
			@param file_names: test the existence of files within the directory
			@type file_names: list[basestring]
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: bool
			@rtype: bool
		"""
		# TODO: test for valid characters

		assert isinstance(silent, bool)
		assert key is None or isinstance(key, basestring)
		assert isinstance(only_parent, bool)
		assert not (only_parent and sub_directories is not None)
		assert not (only_parent and file_names is not None)

		if sub_directories is None:
			sub_directories = []
		if file_names is None:
			file_names = []

		assert directory is None or isinstance(directory, basestring)
		assert isinstance(sub_directories, list)
		assert isinstance(file_names, list)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if directory is None:
			if not silent:
				self._logger.error("{}Invalid directory".format(prefix))
			return False

		if directory == '':
			if not silent:
				self._logger.error("{}Invalid directory: '{}'".format(prefix, directory))
			return False

		directory = self.get_full_path(directory)
		parent_directory = os.path.dirname(directory)
		if not os.path.isdir(parent_directory):
			if not silent:
				self._logger.error("{}Directory does not exist: '{}'".format(prefix, parent_directory))
			return False

		if not only_parent and not os.path.isdir(directory):
			if not silent:
				self._logger.error("{}Directory does not exist: '{}'".format(prefix, directory))
			return False

		for sub_directory in sub_directories:
			if not os.path.isabs(sub_directory):
				sub_directory = os.path.join(directory, sub_directory)
			if not self.validate_dir(sub_directory, key=key):
				return False

		for file_name in file_names:
			if not os.path.isabs(file_name):
				file_name = os.path.join(directory, file_name)
			if not self.validate_file(file_name, key=key):
				return False
		return True

	@staticmethod
	def get_full_path(value):
		"""
			Get the normalized absolute path.

			@attention:

			@param value: directory path or file path
			@type value: basestring

			@return: full path
			@rtype: str
		"""
		assert isinstance(value, basestring)
		value = os.path.expanduser(value)
		value = os.path.normpath(value)
		value = os.path.abspath(value)
		return value

	def validate_number(self, digit, minimum=None, maximum=None, zero=True, key=None, silent=False):
		"""
			Validate that a variable is a number within a specific range if given.

			@attention: valid minimum <= digit <= maximum

			@param digit: Any number such as int, float, long
			@type digit: Number
			@param minimum: valid minimum <= digit
			@type minimum: Number
			@param maximum: valid digit <= maximum
			@type maximum: Number
			@param zero: If 0 is to be excluded
			@type zero: bool
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: bool
			@rtype: bool
		"""
		# TODO: digit >= -1 can not be properly tested yet
		assert isinstance(digit, Number), type(digit)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if minimum and digit < minimum:
			if not silent:
				self._logger.error("{}Invalid digit, must be bigger than {}, but was {}".format(prefix, minimum, digit))
			return False

		if maximum and digit > maximum:
			if not silent:
				self._logger.error("{}Invalid digit, must be smaller than {}, but was {}".format(prefix, maximum, digit))
			return False

		if not zero and digit == 0:
			if not silent:
				self._logger.error("{}Invalid digit, must not be {}".format(prefix, maximum, digit))
			return False
		return True

	def validate_free_space(
		self, directory,
		required_space_in_bytes=None,
		required_space_in_kb=None,
		required_space_in_mb=None,
		required_space_in_gb=None,
		key=None, silent=False
		):
		"""
			Validate that sufficient free space is available at a target directory.

			@attention: Only one 'required space' argument will be accepted

			@param directory: directory path of a folder
			@type directory: basestring
			@param required_space_in_bytes: Required available space in bytes
			@type required_space_in_bytes: Number
			@param required_space_in_kb: Required available space in kilobytes
			@type required_space_in_kb: Number
			@param required_space_in_mb: Required available space in megabytes
			@type required_space_in_mb: Number
			@param required_space_in_gb: Required available space in gigabytes
			@type required_space_in_gb: Number

			@param silent: If True, no error message will be made
			@type silent: bool

			@return: bool
			@rtype: bool
		"""
		required_space = None
		count = 4
		for argument in [required_space_in_bytes, required_space_in_kb, required_space_in_mb, required_space_in_gb]:
			if argument is None:
				count -= 1
			else:
				required_space = argument
		assert count == 1

		# required_space = required_space_in_bytes or required_space_in_kb or required_space_in_mb or required_space_in_gb
		# print required_space, required_space_in_bytes, required_space_in_kb, required_space_in_mb, required_space_in_gb
		assert self.validate_number(required_space, minimum=0)
		assert self.validate_dir(directory, key=key, silent=silent)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		size_label = ""
		free_space = 0
		if required_space_in_bytes is not None:
			size_label = "bytes"
			free_space = self.free_space_in_bytes(directory)

		if required_space_in_kb is not None:
			size_label = "kb"
			free_space = self.free_space_in_kilo_bytes(directory)

		if required_space_in_mb is not None:
			size_label = "mb"
			free_space = self.free_space_in_mega_bytes(directory)

		if required_space_in_gb is not None:
			size_label = "gb"
			free_space = self.free_space_in_giga_bytes(directory)

		if not required_space < free_space:
			if not silent:
				self._logger.error("{}Insufficient space! {:.2f}{label} of {:.2f}{label} available at '{dir}'".format(
					prefix, free_space, required_space, label=size_label, dir=directory))
			return False
		return True

	def free_space_in_giga_bytes(self, directory):
		"""
			Get available free space at a target directory.

			@param directory: directory path of a folder
			@type directory: basestring

			@return: Available free space
			@rtype: float
		"""
		assert self.validate_dir(directory)
		return self._free_space(directory, 3)

	def free_space_in_mega_bytes(self, directory):
		"""
			Get available free space at a target directory.

			@param directory: directory path of a folder
			@type directory: basestring

			@return: Available free space
			@rtype: float
		"""
		assert self.validate_dir(directory)
		return self._free_space(directory, 2)

	def free_space_in_kilo_bytes(self, directory):
		"""
			Get available free space at a target directory.

			@param directory: directory path of a folder
			@type directory: basestring

			@return: Available free space
			@rtype: float
		"""
		assert self.validate_dir(directory)
		return self._free_space(directory, 1)

	def free_space_in_bytes(self, directory):
		"""
			Get available free space at a target directory.

			@param directory: directory path of a folder
			@type directory: basestring

			@return: Available free space
			@rtype: float
		"""
		assert self.validate_dir(directory)
		return self._free_space(directory)

	def _free_space(self, directory, power=0):
		"""
			Get available free space at a target directory.

			@param directory: directory path of a folder
			@type directory: basestring

			@return: Available free space
			@rtype: float
		"""
		assert power >= 0
		assert isinstance(directory, basestring)
		assert self.validate_dir(directory)
		if not directory or not os.path.isdir(directory):
			return 0
		statvfs = os.statvfs(directory)
		free_space = statvfs.f_frsize * statvfs.f_bfree
		return free_space / math.pow(1024, power)
