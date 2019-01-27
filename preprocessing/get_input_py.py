#!/usr/bin/python

'''
	Usege:
	python3 get_input_py.py --input_fasta='./raw_data/PARS_yeast/sce_genes.fasta' 
							--input_score='./raw_data/PARS_yeast/sce_Score.tab' 
							--window_size=37 
							--dataset='py'

	input: sequences with score
	output: top5 and last5 fragments(encode) from each sequence
'''

import os
import argparse
import pandas as pd
import csv

def read_f_fasta(dat_dir):
	with open(dat_dir) as file_obj:
		dat_list = file_obj.readlines()
		dat_list2 = []
		for lines in dat_list:
			lines = lines.strip("\n")
			dat_list2.append(lines)
	return dat_list2

def read_f_score(dat_dir):
	with open(dat_dir) as file_obj:
		dat_list = file_obj.readlines()
		dat_list2 = []
		for lines in dat_list:
			lines = lines.split()
			dat_list2.append(lines)
	return dat_list2


def get_fasta_name(list):
	name_list =[]
	for i in range(len(list)):
		if i%2 == 0:
			name_list.append(list[i][1:])
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

def find_topK_index(your_list,k=3,top='max'):
	if k > len(your_list)-1:
		print("k is larger than length of list")
	else:
		import numpy as np
		arr = np.array(your_list)
		topK_index=[]
		index_seq = np.argsort(arr)

		if top=='max':
			i=0
			while(i>=0):
				if arr[index_seq[len(arr)-k-(i+1)]] ==arr[index_seq[len(arr)-k-i]]:
					i +=1
				else:
					topK_index = index_seq[len(arr)-k-i:]
					break

		if top=='min':
			i=0
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
		#initialize
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

	parser.add_argument('--input_fasta', type=str, default='./raw_data/PARS_yeast/sce_genes.fasta', help='FASTA format file, containing RNA sequences.')
	parser.add_argument('--input_score', type=str, default='./raw_data/PARS_yeast/sce_Score.tab', help='The experimental scores of RNA structural profile.')
	parser.add_argument('--window_size', type=int, default=37, help='The window size when truncating RNA sequences.')
	parser.add_argument('--dataset', type=str, default='py', help='The type of dataset: PARS_human:ph, PARS_yeast:py or NMR_X-Ray:pdb')

	config = parser.parse_args()

	fasta = read_f_fasta(config.input_fasta)#'./sce_genes.fasta'
	score = read_f_score(config.input_score)#'./sce_Score.tab'
	window_len = config.window_size
	
	output_dir = create_dir_if_not_exists(config)
	save_nocode_p = output_dir+'/'+config.dataset+'_nocode_{}.csv'.format(window_len)#./py_nocode_{window_size}
	save_encode_p = output_dir+'/'+config.dataset+'_encode_{}.csv'.format(window_len)#./py_encode_{window_size}
	seq_score = []
	fasta_name_ls=get_fasta_name(fasta)
	score_name_ls=get_score_name(score)

	failed_ls = []
	all_transcript_info = []
	title_lst = ['gene_name','whole_seq','whole_score','pos_index','pos_nucle','pos_frag','neg_index','neg_nucle','neg_frag']
	all_transcript_info.append(title_lst)

	for i in range(len(fasta_name_ls)):
		single_transcript_info=[]
		gene_name = fasta_name_ls[i]
		index_in_score = search_index(gene_name,score_name_ls)
		if index_in_score == None:
			failed_ls.append(gene_name)
		else:
			single_transcript_info.append(gene_name)
			#append() takes exactly one argument
			all_nucle=list(fasta[(2*i)+1]) 
			single_transcript_info.append(all_nucle)
			all_score = list(map(float,list(score[index_in_score][2].split(";"))))
			single_transcript_info.append(all_score)

			max5_index=find_topK_index(all_score,k=5,top='max')
			max5_nucle=[]
			for i in max5_index:
				max5_nucle.append(all_nucle[i])
			max5_seq = get_seq_lst(max5_index,window_len,all_nucle)
			single_transcript_info.append(max5_index)
			single_transcript_info.append(max5_nucle)
			single_transcript_info.append(max5_seq)

			min5_index=find_topK_index(all_score,k=5,top='min')
			min5_nucle=[]
			for i in min5_index:
				min5_nucle.append(all_nucle[i])
			min5_seq = get_seq_lst(min5_index,window_len,all_nucle)
			single_transcript_info.append(min5_index)
			single_transcript_info.append(min5_nucle)
			single_transcript_info.append(min5_seq)

			all_transcript_info.append(single_transcript_info)

	input_data=[]
	feature_lst=['gene_name']+['f'+str(i+1) for i in range(window_len)]+['label','position']
	input_data.append(feature_lst)
	for i in range(1,len(all_transcript_info)):
		for j1 in range(len(all_transcript_info[i][3])):
			single_input=[]
			single_input.append(all_transcript_info[i][0])
			for j1j in range(window_len):
				single_input.append(all_transcript_info[i][5][j1][j1j])
			single_input.append(1)
			single_input.append(all_transcript_info[i][3][j1])
			input_data.append(single_input)
		for j2 in range(len(all_transcript_info[i][6])):
			single_input=[]
			single_input.append(all_transcript_info[i][0])
			for j2j in range(window_len):
				single_input.append(all_transcript_info[i][8][j2][j2j])
			single_input.append(0)
			single_input.append(all_transcript_info[i][6][j2])
			input_data.append(single_input)
	pd.DataFrame(input_data).to_csv(save_nocode_p,header=0,index=0)

	# return all_transcript_info,failed_ls

	#encode
	nocode_data = pd.read_csv(save_nocode_p)
	tocode_data=nocode_data.drop(['gene_name','position'],axis=1)
	nucle_code = {'A':[1,0,0,0],'U':[0,1,0,0],'T':[0,1,0,0],'C':[0,0,1,0],'G':[0,0,0,1],'N':[0,0,0,0],'P':[0,0,0,0],'a':[1,0,0,0],'u':[0,1,0,0],'t':[0,1,0,0],'c':[0,0,1,0],'g':[0,0,0,1],'I':[0,0,0,0]}
	frag_code_input =[]
	code_title = []
	for i in range(1,window_len+2):
		if i<(window_len+1):
			single_title = ['f'+str(i)+'_1','f'+str(i)+'_2','f'+str(i)+'_3','f'+str(i)+'_4']
			code_title.extend(single_title)
		else:
			code_title.append('label')
	

	for i in range(len(tocode_data)):
		one_frag = tocode_data.iloc[i]
		single_code=[]
		for j in range(len(one_frag)-1):
			single_code.extend(nucle_code[one_frag[j]])
		single_code.append(one_frag[-1])
		frag_code_input.append(single_code)

	frag_code_input = pd.DataFrame(frag_code_input)
	gene_name = nocode_data['gene_name']
	full_sequences = (tocode_data.drop('label',axis=1)).apply(lambda x: x.sum(), axis=1)
	frag_code_input.insert(0, 'sequences', full_sequences)
	code_title.insert(0, 'sequences')
	frag_code_input.columns = code_title
	frag_code_input.index = gene_name.tolist()
	frag_code_input.to_csv(save_encode_p,index=True)