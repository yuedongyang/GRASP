#!/usr/bin/python

'''
Usage:
	python3 test.py --input_fasta='test.fa'
					--model_file='./model/py_37.model'
					--window_size=37
					--output_file='./output/predictions'
'''

import argparse
from sklearn.externals import joblib

from utils import *

import datetime
now_time = datetime.datetime.now()
now_time = now_time.strftime("%Y-%m-%d")

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--input_fasta', type=str, required=True, help='Input the test file of RNA sequences(FASTA format, please ensure that your sequences only contain A, C, G, T, and U)')
	parser.add_argument('--model_file', type=str, required=True, help='The full path of the trained model.')
	parser.add_argument('--window_size', type=int, required=True, help='The window size when truncating RNA sequences.')
	parser.add_argument('--output_file', type=str, required=True, help='The full path of the output file.')

	config = parser.parse_args()

	window_size = config.window_size
	input_dict = read_fasta(config.input_fasta)
	check_for_input_sequences(input_dict)

	model_file = config.model_file
	output_file = config.output_file

	predicted_dict = get_features_and_predict(input_dict, window_size, model_file)
	save_file(predicted_dict, input_dict, output_file)