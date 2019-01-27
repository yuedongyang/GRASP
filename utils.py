
import os

import json
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split

import warnings
warnings.filterwarnings("ignore")

dataset_type = {'py':'PARS_yeast', 'ph':'PARS_human', 'pdb':'SS_PDB', 'global':'Consensus'}

def get_params():
	params={
		'booster':['gbtree'],
		'objective':['binary:logistic'],
		'gamma':[0],
		'max_depth':[3,6,9,12],
		'reg_lambda':[0.05,0.1,0.5],
		'reg_alpha':[0.05,0.1,0.5],
		'subsample':[1,0.7],
		'colsample_bytree':[1,0.7],
		'min_child_weight':[1,3,6],
		'silent':[1], 
		'learning_rate':[0.1,0.05],
		'seed':[1000],
		'n_estimators':[500,1000,2000]
	}
	return params

def get_data(file_path):
	with open(file_path, 'r') as f:
		data = pd.read_csv(f, index_col=0)
	return data

def save_csv_data(test_data, test_data_path):
	test_data = pd.DataFrame(test_data)
	test_data.to_csv(test_data_path, index=True)

def get_consensus_data(input_pars_human_file, input_pars_yeast_file, input_pdb_file, config):
	ph_data = get_data(input_pars_human_file)
	py_data = get_data(input_pars_yeast_file)
	pdb_data = get_data(input_pdb_file)

	ph_train, ph_test = train_test_split(ph_data, test_size=0.1, random_state=0)
	py_train, py_test = train_test_split(py_data, test_size=0.1, random_state=0)
	pdb_train, pdb_test = train_test_split(pdb_data, test_size=0.1, random_state=0)

	data_train = ph_train.append(py_train)
	data_train = data_train.append(pdb_train)

	save_test_dir = './data/10per_test'
	if not os.path.exists(save_test_dir):
		os.makedirs(save_test_dir)
	save_csv_data(ph_test, save_test_dir+'/ph_encode_10per_test_{}.csv'.format(config.window_size))
	save_csv_data(py_test, save_test_dir+'/py_encode_10per_test_{}.csv'.format(config.window_size))
	save_csv_data(pdb_test, save_test_dir+'/pdb_encode_10per_test_{}.csv'.format(config.window_size))

	# save_csv_data(ph_train, save_test_dir+'/ph_encode_10per_train_{}.csv'.format(config.window_size))
	# save_csv_data(py_train, save_test_dir+'/py_encode_10per_train_{}.csv'.format(config.window_size))
	# save_csv_data(pdb_train, save_test_dir+'/pdb_encode_10per_train_{}.csv'.format(config.window_size))	

	return data_train, ph_test, py_test, pdb_test

# def get_trained_params(paramsfile):
# 	with open(paramsfile, 'r') as f:
# 		params = json.load(f)
# 	del params['best_score']
# 	return params

def eval(label, preds):
	auc = metrics.roc_auc_score(label, preds)
	acc = metrics.accuracy_score(label, np.round(preds))
	mcc = metrics.matthews_corrcoef(label, np.round(preds))
	ppv = metrics.precision_score(label, np.round(preds))
	recall = metrics.recall_score(label, np.round(preds))
	f1 = metrics.f1_score(label, np.round(preds))
	return auc, acc, mcc, ppv, recall, f1

def create_dir_if_not_exists(dir_path):
	if not os.path.exists(dir_path):
		os.makedirs(dir_path)

def delete_file_if_already_exists(filename):
	if os.path.exists(filename):
		os.remove(filename)

# test.py
def check_for_input_sequences(sequences_dict):
	normarized = ['A', 'G', 'C', 'T', 'U', 'N']
	for seq_name in sequences_dict.keys():
		seq = sequences_dict[seq_name].upper()
		flag = True
		for s in seq:
			if s not in normarized:
				flag = False
				print(s)
				break
		if not flag:
			raise Exception("Invalid RNA Sequences! Please ensure that your sequences only contain A, C, G, T, and U.")

def read_fasta(filename):
	fread = open(filename)
	test_dict = {}
	for line in fread:
		if line.startswith('>'):
			name = line.replace('>', '').split()[0]
			test_dict[name] = ''
		else:
			test_dict[name] += line.replace('\n','')
	fread.close()
	return test_dict

def OneHotEncode(frag_test_dict):
	nucle_code = {'A':[1,0,0,0],'U':[0,1,0,0],'T':[0,1,0,0],'C':[0,0,1,0],'G':[0,0,0,1],'N':[0,0,0,0]}
	# nucle_code = {'A':[1,0,0,0],'U':[0,1,0,0],'T':[0,1,0,0],'C':[0,0,1,0],'G':[0,0,0,1],'N':[0,0,0,0],'P':[0,0,0,0],'a':[1,0,0,0],'u':[0,1,0,0],'t':[0,1,0,0],'c':[0,0,1,0],'g':[0,0,0,1],'I':[0,0,0,0]}
	features = []
	for test_dicts in frag_test_dict:
		fea = []
		for i in test_dicts:
			fea.extend(nucle_code[i])
		features.append(fea)
	return features

def pred_proba(model_name, features):
	clf = joblib.load(model_name)
	preds = clf.predict_proba(features)[:, 1]
	return preds

def get_features_and_predict(input_dict, windows_size, model_name):
	predicted_dict = {}
	# all_features = []

	features_name = []
	for i in range(1, windows_size + 1):
		single_name = ['f'+str(i)+'_1','f'+str(i)+'_2','f'+str(i)+'_3','f'+str(i)+'_4']
		features_name.extend(single_name)

	for key in input_dict.keys():
		# print('{} is processing'.format(key))
		raw_test_dict = input_dict[key]
		slen = len(raw_test_dict)

		all_frag_test_dict = []
		for idx in range(slen):

			#get sub-test_dictuence
			side_len = int((windows_size - 1) / 2)
			init_index = idx - side_len
			end_index = idx + side_len
			frag_test_dict = []
			for i in range(windows_size):
				frag_test_dict.append('N')

			if (init_index < 0 and end_index <= slen - 1):
				for i in range(end_index + 1):
					frag_test_dict[abs(init_index) + i] = raw_test_dict[i]
			elif (init_index >= 0 and end_index > slen - 1):
				for i in range(slen - init_index):
					frag_test_dict[i] = raw_test_dict[init_index + i]
			elif (init_index < 0 and end_index > slen - 1):
				for i in range(slen):
					frag_test_dict[abs(init_index) + i] = raw_test_dict[i]
			else:
				for i in range(windows_size):
					frag_test_dict[i] = raw_test_dict[init_index + i]

			all_frag_test_dict.append(frag_test_dict)

		features = OneHotEncode(all_frag_test_dict)

		features = pd.DataFrame(features, columns=features_name)
		if key == 'seq1':
			print(features)

		#predict
		preds = pred_proba(model_name, features)
		predicted_dict[key] = preds

	# all_features = pd.DataFrame(all_features)
	# all_features.to_csv('test_features.csv', index=False)
	return predicted_dict

def save_file(predicted_dict, input_dict, filename):
	fwrite = open(filename, 'w')
	for key in input_dict.keys():
		fwrite.write('>')
		fwrite.write(key)
		fwrite.write('\n')
		fwrite.write(input_dict[key])
		fwrite.write('\n')
		fwrite.write(";".join(str(i) for i in predicted_dict[key]))
		fwrite.write(";")
		fwrite.write('\n')
	fwrite.close()
	print("Saving predictions: ", filename)