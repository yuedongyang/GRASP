#!/usr/bin/python

'''
	Usege:
	python3 get_input_ph.py --input_fasta='./raw_data/NMR_X-Ray/new75.fa' 
							--input_ss='./raw_data/NMR_X-Ray/pdball.SS' 
							--window_size=37 
							--dataset='pdb'

	input: sequences with score
	output: top5 and last5 fragments(encode) from each sequence
'''
import os
import argparse
import pandas as pd
import csv


def rd_normal(dat_dir):
	with open(dat_dir) as file_obj:
		dat_list = file_obj.readlines()
	return dat_list

def rd_fasta(dat_dir):
	with open(dat_dir) as file_obj:
		dat_list = file_obj.readlines()
		dat_list2 = []
		for i in range(0,len(dat_list),2):
			dat_list2.append([dat_list[i],dat_list[i+1]])
	return dat_list2

def rd_score(dat_dir):
	with open(dat_dir) as file_obj:
		#dat_list = file_obj.readlines()
		dat_list2 = []
		for lines in file_obj:
			dat_list2.append(lines.split())
	return dat_list2

def rd_fasta_name(dat_dir):
	with open(dat_dir) as file_obj:
		dat_list = file_obj.readlines()
		dat_list2 = []
		for i in range(0,len(dat_list),2):
			dat_list2.append(dat_list[i][1:-1])
	return dat_list2

def get_fasta_name(list):
	name_list =[]
	for i in range(0,len(list)):
		name_list.append(list[i][0][1:-2])
	return name_list

def get_score_name(list):
	name_list =[]
	for i in range(len(list)):
		name_list.append(list[i][0])
	return name_list

def search_index(gene_name,list):
	try:
		gene_index = list.index(gene_name)
	except ValueError:
		return None
	else:
		return gene_index

def find_topK_index(your_list,k=5,top='max'):
	if k > len(your_list)-1:
		print("k is larger than length of list")
	else:
		import numpy as np
		arr = np.array(your_list)
		topK_index=[]
		index_seq = np.argsort(arr)

		if top=='max':
			if abs(arr[index_seq[-1]])==0.0:
				return None
			else:
				i=0
				if abs(arr[index_seq[len(arr)-k]])== 0.0:
					while abs(arr[index_seq[len(arr)-k+i]])== 0.0:
						i +=1
					topK_index = index_seq[len(arr)-k+i:]
				else:
					while(i>=0):
						if arr[index_seq[len(arr)-k-(i+1)]] ==arr[index_seq[len(arr)-k-i]]:
							i +=1
						else:
							topK_index = index_seq[len(arr)-k-i:]
							break

		if top=='min':
			if abs(arr[index_seq[0]])==0.0:
				return None
			else:
				i=0
				if abs(arr[index_seq[k-1]])== 0.0:
					while abs(arr[index_seq[k-1-i]])== 0.0:
						i +=1
					topK_index = index_seq[:(k-1-i+1)]
				else:
					while(i>=0):
						if arr[index_seq[k+i]] ==arr[index_seq[(k-1)+i]]:
							i +=1
						else:
							topK_index = index_seq[:(k+i)]
							break
		return list(topK_index)

def find_pair_index(your_list,top='pair'):
	pair_index=list()
	#index_seq = np.argsort(arr)
	if top=='pair':
		pair_index = [i for i in range(len(your_list)) if your_list[i] >0.0]
	if top=='unpair':
		pair_index = [i for i in range(len(your_list)) if your_list[i] == 0.0]
	return list(pair_index)

def get_seq_lst(cenucle_idx_lst,frag_len,seq_lst):
	side_len = int((frag_len-1)/2)
	all_frag_seq = []
	for idx in cenucle_idx_lst:
		init_idx = idx-side_len
		end_idx = idx+side_len
		frag_seq = []
		for i in range(frag_len):
			frag_seq.append('N')

		if init_idx<0 and end_idx<=(len(seq_lst)-1):
			for i in range(end_idx+1):
				frag_seq[abs(init_idx)+i]=seq_lst[i]
		elif init_idx>=0 and end_idx>(len(seq_lst)-1):
			for i in range(len(seq_lst)-init_idx):
				frag_seq[i]=seq_lst[init_idx+i]
		elif init_idx<0 and end_idx>(len(seq_lst)-1):
			for i in range(len(seq_lst)):
				frag_seq[abs(init_idx)+i]=seq_lst[i]
		else:
			for i in range(frag_len):
				frag_seq[i]=seq_lst[init_idx+i]
		all_frag_seq.append(frag_seq)
	return all_frag_seq

dataset_type = {'py':'PARS_yeast', 'ph':'PARS_human', 'pdb':'NMR_X-Ray'}
def create_dir_if_not_exists(config):
	dir_path = '../data/' + dataset_type[config.dataset]
	if not os.path.exists(dir_path):
		os.makedirs(dir_path)
	return dir_path

if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument('--input_fasta', type=str, default='./raw_data/NMR_X-Ray/new75.fa', help='FASTA format file, containing RNA sequences.')
	parser.add_argument('--input_ss', type=str, default='./raw_data/NMR_X-Ray/pdball.SS', help='The RNA sequences structural profile.')
	parser.add_argument('--window_size', type=int, default=37, help='The window size when truncating RNA sequences.')
	parser.add_argument('--dataset', type=str, default='pdb', help='The type of dataset: PARS_human:ph, PARS_yeast:py or NMR_X-Ray:pdb')

	config = parser.parse_args()


	pdb_select_p=str(config.input_fasta)	#new75.fa
	pdb_all_p=str(config.input_ss)		#pdball.SS
	window_len=int(config.window_size)		#37
	
	output_dir = create_dir_if_not_exists(config)
	save_nocode_p = output_dir+'/'+config.dataset+'_nocode_{}.csv'.format(window_len)#./py_nocode_{window_size}
	save_encode_p = output_dir+'/'+config.dataset+'_encode_{}.csv'.format(window_len)#./py_encode_{window_size}
	
	select_name = rd_fasta_name(pdb_select_p)
	pdb_all_lst = rd_normal(pdb_all_p)

	all_data_info =list()
	all_data_info.append(['gene_name','whole_seq','whole_score','pos_frag','pos_ix','neg_frag','neg_ix'])
	for i in range(0,len(pdb_all_lst),3):
		tmp_name = pdb_all_lst[i][1:-1]
		if tmp_name in select_name:
			single_data_info = list()
			single_data_info.append(pdb_all_lst[i][1:-1])
			a_seq=pdb_all_lst[i+1]
			a_score=pdb_all_lst[i+2]
			single_data_info.append(a_seq)
			single_data_info.append(a_score)
			pos_seq_ix = find_pair_index(list(map(float, a_score.split())),top='pair')
			pos_seq_lst = get_seq_lst(pos_seq_ix,window_len,seq_lst=list(a_seq[:-1]))
			single_data_info.append(pos_seq_lst)
			single_data_info.append(pos_seq_ix)
			neg_seq_ix = find_pair_index(list(map(float, a_score.split())),top='unpair')
			neg_seq_lst = get_seq_lst(neg_seq_ix,window_len,seq_lst=list(a_seq[:-1]))
			single_data_info.append(neg_seq_lst)
			single_data_info.append(neg_seq_ix)
			all_data_info.append(single_data_info)
	# print('get seq over')

	input_data=[]
	feature_lst=['gene_name']+['f'+str(i+1) for i in range(window_len)]+['label','position']
	input_data.append(feature_lst)
	for i in range(1,len(all_data_info)):
		for j1 in range(len(all_data_info[i][3])):
			single_input=[]
			single_input.append(all_data_info[i][0])
			for j1j in range(window_len):
				single_input.append(all_data_info[i][3][j1][j1j])
			single_input.append(1)
			single_input.append(all_data_info[i][4][j1])
			input_data.append(single_input)
		for j2 in range(len(all_data_info[i][5])):
			single_input=[]
			single_input.append(all_data_info[i][0])
			for j2j in range(window_len):
				single_input.append(all_data_info[i][5][j2][j2j])
			single_input.append(0)
			single_input.append(all_data_info[i][6][j2])
			input_data.append(single_input)
				
	pd.DataFrame(input_data).to_csv(save_nocode_p,header=0,index=0)


	#encode
	PDB_input_data = pd.read_csv(save_nocode_p)
	n_PDB_input_data=PDB_input_data.drop(['gene_name','position'],axis=1)
	nucle_code = {'A':[1,0,0,0],'U':[0,1,0,0],'T':[0,1,0,0],'C':[0,0,1,0],'G':[0,0,0,1],'N':[0,0,0,0],'P':[0,0,0,0],'a':[1,0,0,0],'u':[0,1,0,0],'t':[0,1,0,0],'c':[0,0,1,0],'g':[0,0,0,1],'I':[0,0,0,0]}
	frag_code_input =[]
	code_title = []
	for i in range(1,window_len+2):
		if i<(window_len+1):
			single_title = ['f'+str(i)+'_1','f'+str(i)+'_2','f'+str(i)+'_3','f'+str(i)+'_4']
			code_title.extend(single_title)
		else:
			code_title.append('label')
	
	for i in range(len(n_PDB_input_data)):
		one_frag = n_PDB_input_data.iloc[i]
		single_code=[]
		for j in range(len(one_frag)-1):
			single_code.extend(nucle_code[one_frag[j]])
		single_code.append(one_frag[-1])
		frag_code_input.append(single_code)

	frag_code_input = pd.DataFrame(frag_code_input)
	gene_name = PDB_input_data['gene_name']
	full_sequences = (n_PDB_input_data.drop('label', axis=1)).apply(lambda x: x.sum(), axis=1)
	frag_code_input.insert(0, 'sequences', full_sequences)
	code_title.insert(0, 'sequences')
	frag_code_input.columns = code_title
	frag_code_input.index = gene_name.tolist()
	frag_code_input.to_csv(save_encode_p,index=True)