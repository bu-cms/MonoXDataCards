from hepdata_validator.data_file_validator import DataFileValidator
import argparse

parser = argparse.ArgumentParser(description='Validate yaml files.')
parser.add_argument('-filename',dest='filename', type=str, help='file to check')

args = parser.parse_args()

data_file_validator = DataFileValidator()

# the validate method takes a string representing the file path.
data_file_validator.validate(file_path=args.filename)

# if there are any error messages, they are retrievable through this call
data_file_validator.get_messages()

# the error messages can be printed
data_file_validator.print_errors('data.yaml')
