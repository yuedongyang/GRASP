#!/usr/bin/python

'''
Usage:
	python3 evaluation.py --mode='single'
						  --test_dataset='ph'
						  --train_model='py'
						  --window_size=37
'''

import argparse
import json
from sklearn.externals import joblib

from utils import *


if __name__=="__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument('--mode', type=str, default='single', help='The trained model is trained by two training modes: single or global, \'single\' is trained by single dataset and \'global\' is trained by all datasets.')
	parser.add_argument('--test_dataset', type=str, default='ph', help='This is for test-dataset. Choose from: \'py\' is PARS-Yeast \'ph\' is PARS-Human, \'pdb\' is NMR/X-ray.')
	parser.add_argument('--train_model', type=str, default='py', help='This is for using which trained model. Choose from: \'py\' is PARS-Yeast \'ph\' is PARS-Human, \'pdb\' is NMR/X-ray and \'global\' is Consensus-model.')
	parser.add_argument('--window_size', type=int, default=37, help='The window size when truncating RNA sequences.')
	parser.add_argument('--output', type=str, default='./preds/', help='The directory for saving predictions on test-dataset.')

	config = parser.parse_args()	

	window_size = config.window_size
	outputfile_dir = config.output + dataset_type[config.test_dataset]
	create_dir_if_not_exists(outputfile_dir)

	if (config.mode == 'single'):
		inputfile = './data/' + dataset_type[config.test_dataset] + '/' + config.test_dataset + '_encode_{}.csv'.format(window_size)
		model = './model/' + dataset_type[config.train_model] + '/' + config.train_model +'_{}.model'.format(window_size)
		outputfile = outputfile_dir + '/' + config.test_dataset + '_predicted_by_' + config.train_model + '_{}.csv'.format(window_size)
		
		input_data = get_data(inputfile)
		output = input_data
		label_test = input_data['label']
		input_data = input_data.drop(['label', 'sequences'], axis=1)

		clf = joblib.load(model)
		preds = clf.predict_proba(input_data)[:,1]

		output['preds'] = preds
		# output.to_csv(outputfile, index=True)

		print ('training model:', config.train_model, 'test dataset:', config.test_dataset)
		auc, acc, mcc, ppv, recall, f1 = eval(label_test, preds)
		print ('auc:', auc)
		print ('acc:', acc)
		print ('mcc:', mcc)
		print ('ppv:', ppv)
		print ('recall:', recall)
		print ('f1-score:', f1)

	else :
		input_pars_human_file = './data/10per_test/ph_encode_10per_test_{}.csv'.format(window_size)
		input_pars_yeast_file = './data/10per_test/py_encode_10per_test_{}.csv'.format(window_size)
		input_pdb_file = './data/10per_test/pdb_encode_10per_test_{}.csv'.format(window_size)

		output_pars_human_preds_file = outputfile_dir + '/ph_encode_10per_test_preds_{}.csv'.format(window_size)
		output_pars_yeast_preds_file = outputfile_dir + '/py_encode_10per_test_preds_{}.csv'.format(window_size)
		output_pdb_preds_file = outputfile_dir + '/pdb_encode_10per_test_preds_{}.csv'.format(window_size)

		model = './model/' + dataset_type[config.train_model] + '/' + config.train_model +'_{}.model'.format(window_size)

		pars_human_test_data = get_data(input_pars_human_file)
		pars_yeast_test_data = get_data(input_pars_yeast_file)
		pdb_test_data = get_data(input_pdb_file)

		ph_test_preds = pars_human_test_data
		py_test_preds = pars_yeast_test_data
		pdb_test_preds = pdb_test_data

		pars_human_test_label = pars_human_test_data['label']
		pars_yeast_test_label = pars_human_test_data['label']
		pdb_test_label = pdb_test_data['label']

		pars_human_test_data = pars_human_test_data.drop(['label', 'sequences'], axis=1)
		pars_human_test_data = pars_yeast_test_data.drop(['label', 'sequences'], axis=1)
		pdb_test_data = pdb_test_data.drop(['label', 'sequences'], axis=1)

		clf = joblib.load(model)
		ph_preds = clf.predict_proba(pars_human_test_data)[:,1]
		py_preds = clf.predict_proba(pars_yeast_test_data)[:,1]
		pdb_preds = clf.predict_proba(pdb_test_data)[:,1]

		ph_test_preds['preds'] = ph_preds
		ph_test_preds.to_csv(output_pars_human_preds_file, index=True)
		py_test_preds['preds'] = py_preds
		py_test_preds.to_csv(output_pars_yeast_preds_file, index=True)
		pdb_test_preds['preds'] = pdb_preds
		pdb_test_preds.to_csv(output_pdb_preds_file, index=True)

		print ('training model:', config.train_model)
		py_auc = metrics.roc_auc_score(pars_yeast_test_label, py_preds)
		print ('auc of test dataset: PARS-Yeast:', py_auc)
		ph_auc = metrics.roc_auc_score(pars_human_test_label, pars_human_test_label)
		print ('auc of test dataset: PARS-Human:', ph_auc)
		pdb_auc = metrics.roc_auc_score(pdb_test_label, pdb_preds)
		print ('auc of test dataset: NMR/X-ray:', pdb_auc)

		# save_auc = 'score.json'
		# score_json = {}
		# score_json['ph'] = roc_auc_ph
		# score_json['py'] = roc_auc_py
		# score_json['pdb'] = roc_auc_pdb
		# with open(save_auc,'a') as f:
		# 	json.dump(score_json, f, ensure_ascii=False)
		# 	f.write('\n')
		# print("save end")
	