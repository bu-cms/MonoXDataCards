from hepdata_validator.submission_file_validator import SubmissionFileValidator

submission_file_validator = SubmissionFileValidator()
submission_file_path = 'submission/submission.yaml'

# the validate method takes a string representing the file path.
is_valid_submission_file = submission_file_validator.validate(file_path=submission_file_path)

# if there are any error messages, they are retrievable through this call
submission_file_validator.get_messages()

# the error messages can be printed
submission_file_validator.print_errors(submission_file_path)
