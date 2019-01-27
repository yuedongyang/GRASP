#!/usr/bin/python


'''
Usage:
	python3 train.py --input='./data/PARS_yeast/py_encode_37.csv'
					 --training_mode='single'
					 --dataset='py'
					 --window_size=37
					 --n_jobs=10

'''

import json
import argparse
import xgboost as xgb
from sklearn.externals import joblib
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV

from utils import *

import warnings
warnings.filterwarnings("ignore")

if __name__=="__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument('--input', type=str, default='./data/PARS_human/ph_encode_37.csv', help='The full path of the input file.')
	parser.add_argument('--training_mode', type=str, default='single', help='There are two training mode: single or global, \'single\' is trained by single dataset and \'global\' is trained by all datasets.')
	parser.add_argument('--dataset', type=str, default='ph', help='This is for dataset, \'py\' is for PARS-Yeast, \'ph\' is for PARS-Human, \'pdb\' is for NMR/X-ray and \'global\' is for the three dataset mixed.')
	parser.add_argument('--window_size', type=int, default=37, help='The window size when truncating RNA sequences.')
	parser.add_argument('--n_jobs', type=int, default=-1, help='Number of jobs to run in parallel.')

	config = parser.parse_args()

	window_size = config.window_size

	paras_dir_path = './parameters/' + dataset_type[config.dataset]
	model_dir_path = './model/' + dataset_type[config.dataset]
	create_dir_if_not_exists(paras_dir_path)
	create_dir_if_not_exists(model_dir_path)
	save_best_parameters = paras_dir_path + '/' + config.dataset +'_{}.json'.format(window_size)
	save_best_model = model_dir_path + '/' + config.dataset +'_{}.model'.format(window_size)

	if (config.training_mode == 'single'):
		if (config.dataset not in ['py', 'ph', 'pdb']):
			raise Exception("Invalid dataset's name for single training!", config.dataset, "Please choose from: \'py\', \'ph\' and \'pdb\'.")
		# inputfile='./data/' + dataset_type[config.dataset] + '/' + config.dataset + '_encode_{}.csv'.format(window_size)
		inputfile = config.input

		train_data = get_data(inputfile)

	elif (config.training_mode == 'global'):
		input_pars_human_file='./data/PARS_human/ph_encode_{}'.format(window_size)+'.csv'
		input_pars_yeast_file='./data/PARS_yeast/py_encode_{}'.format(window_size)+'.csv'
		input_pdb_file='./data/NMR_X-Ray/pdb_encode_{}'.format(window_size)+'.csv'

		train_data, ph_test, py_test, pdb_test = get_consensus_data(input_pars_human_file, input_pars_yeast_file, input_pdb_file, config)
		
	else:
		raise Exception("Invalid Training mode! Please choose from \'single\' and \'global\'.")

	label_train = train_data['label']
	train_data = train_data.drop('label', axis=1)
	train_data = train_data.drop('sequences', axis=1)

	xgb_train = xgb.DMatrix(train_data, label=label_train)
	xgb_model = xgb.XGBClassifier()
	parameters = get_params()
	print('Begin training with dataset: ', config.dataset)
	clf = GridSearchCV(xgb_model, param_grid=parameters ,n_jobs=config.n_jobs, cv=StratifiedKFold(label_train,n_folds=5,shuffle=True), scoring='roc_auc', verbose=2, refit=True)
	clf.fit(train_data,label_train)

	best_parameters, score, _ = max(clf.grid_scores_, key=lambda x:x[1])
	print('Raw AUC score:',score)
	for param_name in sorted(best_parameters.keys()):
		print("%s:%r"%(param_name, best_parameters[param_name]))
	best_parameters['best_score'] = score

	delete_file_if_already_exists(save_best_parameters)
	with open(save_best_parameters, 'a') as f:
		json.dump(best_parameters, f, ensure_ascii=False)
		f.write('\n')
	print ("Save best parameters: ", save_best_parameters)

	joblib.dump(clf, save_best_model)
	print ("Save model: ", save_best_model)