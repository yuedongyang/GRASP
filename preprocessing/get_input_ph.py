#!/usr/bin/python


'''
	Usege:
	python3 get_input_ph.py --input_fasta='./raw_data/PARS_human/pars_h_gene_haveseq.fasta' 
							--input_score='./raw_data/PARS_human/GM12878_score_no+5score_haveSeq.tab' 
							--window_size=37 
							--dataset='ph'

	input: sequences with score
	output: top5 and last5 fragments(encode) from each sequence
'''
import os
import argparse
import pandas as pd
import csv


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

def get_fasta_name(list):
	name_list =[]
	for i in range(0,len(list)):
		name_list.append(list[i][0][1:-1])
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

def get_values(ix_lst,lst):
	choose_value=[]
	for ix in ix_lst:
		choose_value.append(lst[ix])
	return choose_value

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

	parser.add_argument('--input_fasta', type=str, default='./raw_data/PARS_human/pars_h_gene_haveseq.fasta', help='FASTA format file, containing RNA sequences.')
	parser.add_argument('--input_score', type=str, default='./raw_data/PARS_human/GM12878_score_no+5score_haveSeq.tab', help='The experimental scores of RNA structural profile.')
	parser.add_argument('--window_size', type=int, default=37, help='The window size when truncating RNA sequences.')
	parser.add_argument('--dataset', type=str, default='ph', help='The type of dataset: PARS_human:ph, PARS_yeast:py or NMR_X-Ray:pdb')

	config = parser.parse_args()


	fasta_p = (config.input_fasta)
	score_p = (config.input_score)
	window_len = config.window_size
	
	output_dir = create_dir_if_not_exists(config)
	save_nocode_p = output_dir+'/'+config.dataset+'_nocode_{}.csv'.format(window_len)#./ph_nocode_{window_size}
	save_encode_p = output_dir+'/'+config.dataset+'_encode_{}.csv'.format(window_len)#./ph_encode_{window_size}

	fasta = rd_fasta(fasta_p);			#print(fasta[0])
	score = rd_score(score_p);			#print(score[0])
	seq_score = []
	fasta_name_ls=get_fasta_name(fasta);#print(fasta_name_ls[:2])
	score_name_ls=get_score_name(score);#print(score_name_ls[:2])
	failed_ls = []

	all_data_info =list()
	all_data_info.append(['gene_name','whole_seq','whole_score','pos_frag','pos_ix','neg_frag','neg_ix'])

	for g_s in score:
		single_data_info = list()
		a_gene_name = g_s[0];#print(a_gene_name)
		a_score = list(map(float,g_s[1].split(';')[:-1]))
		a_seq = list(fasta[fasta_name_ls.index(str(a_gene_name))][1])[:-1]
		single_data_info.append(a_gene_name)
		single_data_info.append(a_seq)
		single_data_info.append(a_score)
		pos_seq_ix = find_topK_index(a_score,top='max')
		#pos_seq_ix_values = get_values(pos_seq_ix,a_score)
		if pos_seq_ix != None:
			pos_seq_ix_values = get_values(pos_seq_ix,a_score)
			pos_seq_lst = get_seq_lst(pos_seq_ix,window_len,seq_lst=a_seq)
			single_data_info.append(pos_seq_lst)
			single_data_info.append(pos_seq_ix)
		else:
			single_data_info.append([])
			single_data_info.append([])
		neg_seq_ix = find_topK_index(a_score,top='min')
		#neg_seq_ix_values = get_values(neg_seq_ix,a_score)
		if neg_seq_ix != None:
			neg_seq_ix_values = get_values(neg_seq_ix,a_score)
			neg_seq_lst = get_seq_lst(neg_seq_ix,window_len,seq_lst=a_seq)
			single_data_info.append(neg_seq_lst)
			single_data_info.append(neg_seq_ix)
		else:
			single_data_info.append([])
			single_data_info.append([])
		all_data_info.append(single_data_info)
	# print('get seq over')

	input_data=[]
	feature_lst=['gene_name']+['f'+str(i+1) for i in range(window_len)]+['label','position']
	input_data.append(feature_lst)
	for i in range(1,len(all_data_info)):
		if all_data_info[i][3] != []:
			for j1 in range(len(all_data_info[i][3])):
				single_input=[]
				single_input.append(all_data_info[i][0])
				for j1j in range(window_len):
					single_input.append(all_data_info[i][3][j1][j1j])
				single_input.append(1)
				single_input.append(all_data_info[i][4][j1])
				input_data.append(single_input)
		if all_data_info[i][5] != []:
			for j2 in range(len(all_data_info[i][5])):
				single_input=[]
				single_input.append(all_data_info[i][0])
				for j2j in range(window_len):
					single_input.append(all_data_info[i][5][j2][j2j])
				single_input.append(0)
				single_input.append(all_data_info[i][-1][j2])
				input_data.append(single_input)
	#similar value choose
	all_frag = list()
	for i in range(1,len(input_data)):
		all_frag.append([''.join(input_data[i][1:-2]),input_data[i][-2]])
	all_frag_df = pd.DataFrame(all_frag)
	new_input_data=[]
	feature_lst=['gene_name']+['f'+str(i+1) for i in range(window_len)]+['label','position']
	new_input_data.append(feature_lst)
	appended_frag=list()
	for i in range(len(all_frag_df)):
		# if i%1000 ==0:
		# 	print(i)
		target_frag = all_frag_df.iloc[i,0]
		if target_frag in appended_frag:
			continue
		target_frag_label = int(all_frag_df.iloc[i,1])
		tmp_df = all_frag_df.loc[all_frag_df[0]==target_frag]
		if len(tmp_df) > 1 :
			tmp_df2 = tmp_df.loc[tmp_df[1]==target_frag_label]
			if len(tmp_df2) >= 0.9*len(tmp_df):
				appended_frag.append(target_frag)
				new_input_data.append(input_data[i+1])

	pd.DataFrame(new_input_data).to_csv(save_nocode_p,header=0,index=0)

	#encode
	PARS_input_data = pd.read_csv(save_nocode_p)
	n_PARS_input_data=PARS_input_data.drop(['gene_name','position'],axis=1)
	nucle_code = {'A':[1,0,0,0],'U':[0,1,0,0],'T':[0,1,0,0],'C':[0,0,1,0],'G':[0,0,0,1],'N':[0,0,0,0],'P':[0,0,0,0],'a':[1,0,0,0],'u':[0,1,0,0],'t':[0,1,0,0],'c':[0,0,1,0],'g':[0,0,0,1],'I':[0,0,0,0]}
	frag_code_input =[]
	code_title = []
	for i in range(1,window_len+2):
		if i<(window_len+1):
			single_title = ['f'+str(i)+'_1','f'+str(i)+'_2','f'+str(i)+'_3','f'+str(i)+'_4']
			code_title.extend(single_title)
		else:
			code_title.append('label')

	for i in range(len(n_PARS_input_data)):
		one_frag = n_PARS_input_data.iloc[i]
		single_code=[]
		for j in range(len(one_frag)-1):
			single_code.extend(nucle_code[one_frag[j]])
		single_code.append(one_frag[-1])
		frag_code_input.append(single_code)

	frag_code_input = pd.DataFrame(frag_code_input)
	gene_name = PARS_input_data['gene_name']
	full_sequences = (n_PARS_input_data.drop('label',axis=1)).apply(lambda x: x.sum(), axis=1)
	frag_code_input.insert(0, 'sequences', full_sequences)
	code_title.insert(0, 'sequences')
	frag_code_input.columns = code_title
	frag_code_input.index = gene_name.tolist()
	frag_code_input.to_csv(save_encode_p,index=True)