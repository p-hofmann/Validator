__author__ = 'hofmann'
__version__ = '0.0.1'

import io
import StringIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from validator import Validator


class SequenceValidator(Validator):
	"""
		Colection of methods for the validation of sequences and sequence files
	"""
	_formats = ["fasta", "fastq"]

	_sequence_indicators = {
		"fasta": ">",
		"fastq": "@"
		}

	_alphabets = {
		"rna": [IUPAC.unambiguous_rna, IUPAC.ambiguous_rna],
		"dna": [IUPAC.unambiguous_dna, IUPAC.ambiguous_dna, IUPAC.extended_dna],
		"protein": [IUPAC.protein, IUPAC.extended_protein]
		}

	_illegal_characters = {
		'\r': "carriage return",
		'\0': "null character"
		}

	@staticmethod
	def _is_stream(stream):
		"""
			Test for streams

			@param stream: Any kind of stream type
			@type stream: file | io.FileIO | StringIO.StringIO

			@return: True if stream
			@rtype: bool
		"""
		return isinstance(stream, (file, io.FileIO, StringIO.StringIO)) or stream.__class__ is StringIO.StringIO

	def validate_sequence_file(self, file_path, file_format, sequence_type, ambiguous, key=None, silent=False):
		"""
			Validate a file to be correctly formatted

			@param file_path: Path to file containing sequences
			@type file_path: str | unicode
			@param file_format: Format of the file at the file_path provided. Valid: 'fasta', 'fastq'
			@type file_format: str | unicode
			@param sequence_type: Are the sequences DNA or RNA? Valid: 'rna', 'dna', 'protein'
			@type sequence_type: str | unicode
			@param ambiguous: True or False, DNA example for strict 'GATC',  ambiguous example 'GATCRYWSMKHBVDN'
			@type ambiguous: bool
			@param key: If True, no error message will be made
			@type key: basestring | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if the file is correctly formatted
			@rtype: bool
		"""
		assert self.validate_file(file_path)
		assert isinstance(file_format, basestring)
		file_format = file_format.lower()
		assert file_format in self._formats
		assert isinstance(sequence_type, basestring)
		sequence_type = sequence_type.lower()
		assert sequence_type in self._alphabets

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if ambiguous:
			alphabet = self._alphabets[sequence_type][1]
		else:
			alphabet = self._alphabets[sequence_type][0]

		with open(file_path) as file_handle:
			if not self._validate_file_start(file_handle, file_format):
				if not silent:
					self._logger.error("{}Invalid begin of file.".format(prefix))
				return False
			for seq_record in SeqIO.parse(file_handle, file_format, alphabet=alphabet):
				result = self.validate_sequence(seq_record.seq, key=key, silent=silent)
				if not result:
					if not silent:
						self._logger.error("{}Invalid sequence: '{}'.".format(prefix, seq_record.id))
					return False
		return True

	def _validate_file_start(self, file_handle, file_format):
		"""
			Validate that a stream with sequences starts with the correct character

			@param file_handle: Any kind of stream type
			@type file_handle: file | io.FileIO | StringIO.StringIO
			@param file_format: Format of the file at the file_path provided. Valid: 'fasta', 'fastq'
			@type file_format: str | unicode

			@return: True if the first character is correct
			@rtype: bool
		"""
		assert self._is_stream(file_handle)
		assert isinstance(file_format, basestring)
		file_format = file_format.lower()
		assert file_format in self._formats

		sequence_indicator = self._sequence_indicators[file_format]

		head = file_handle.read(1)
		file_handle.seek(0)
		if not head:
			return False
		if not head.startswith(sequence_indicator):
			return False
		return True

	def validate_sequence_id(self, identifier, key=None, silent=False):
		"""
			Validate that the sequence identifier has only valid characters

			@attention:

			@param identifier: sequence
			@type identifier: str | unicode
			@param key: If True, no error message will be made
			@type key: str | unicode | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		assert isinstance(identifier, basestring)
		assert isinstance(silent, bool)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if not len(identifier) > 0:
			if not silent:
				self._logger.error("{}Missing sequence id".format(prefix))
			return False
		for illegal_character in self._illegal_characters:
			if illegal_character in identifier:
				if not silent:
					self._logger.error("{}Illegal {} in sequence id".format(
						prefix, self._illegal_characters[illegal_character]))
				return False
		return True

	def validate_sequence_description(self, description, key=None, silent=False):
		"""
			Validate that the sequence description has only valid characters

			@attention:

			@param description: sequence
			@type description: str | unicode
			@param key: If True, no error message will be made
			@type key: str | unicode | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		assert isinstance(description, basestring)
		assert isinstance(silent, bool)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		for illegal_character in self._illegal_characters:
			if illegal_character in description:
				if not silent:
					self._logger.error("{}Illegal {} in sequence description".format(
						prefix, self._illegal_characters[illegal_character]))
				return False
		return True

	def validate_sequence(self, sequence, key=None, silent=False):
		"""
			Validate that the sequence has only valid characters

			@attention:

			@param sequence: sequence
			@type sequence: Seq
			@param key: If True, no error message will be made
			@type key: basestring | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		assert isinstance(sequence, Seq)
		assert isinstance(silent, bool)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		set_alphabet = set(sequence.alphabet.letters)
		if not len(sequence) > 0:
			if not silent:
				self._logger.error("{}Empty sequence".format(prefix))
			return False
		if not set_alphabet.issuperset(sequence.upper()):
			if not silent:
				difference = set(sequence.upper()).difference(set_alphabet)
				difference.discard(set_alphabet)
				self._logger.error("{}Invalid characters: '{}'".format(prefix, ", ".join(difference)))
			return False
		return True