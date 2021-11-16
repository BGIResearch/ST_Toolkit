#!/usr/bin/env python3

##Author: xujunhao xujunhao@genomics.cn
##Date: 2020-6-1
##modified by liuxing on 2021-09-18

import os
from optparse import OptionParser
import sys
from collections import defaultdict
import re
import json
import threading
from numpy import *

def main():
	"""
	%prog [options]
	stat spatial filter result to visualization report
	"""

	parser = OptionParser(main.__doc__)
	parser.add_option("-p", "--path_result", help = "The filter path")
	parser.add_option("-o", "--outfile" , default = None, help = "The output file")
	opts, args = parser.parse_args()

	if(opts.path_result == None):
		sys.exit(not parser.print_help())

	result_path = opts.path_result
	outfile_path = opts.outfile

	if outfile_path == None:
		 outfile_path = result_path + '/'

	Statistic = statistic(result_path, outfile_path)
	Statistic.run()

class statistic():
	def __init__(self, result_path, outfile_path):
		self.result_path = result_path
		self.outfile_path = outfile_path
		self.alignment_result_path = result_path + "/result/"
		self.exp_result_path = result_path + '/result/GetExp/'
		self.filter_result_path = result_path + "/result/"
		self.figure_path = result_path + '/result/GetExp/tissueCut/segmentation/TissueFig/'
		self.cell_cluster = result_path + '/result/GetExp/tissueCut/cell_cluster.h5ad'
		self.tissueCut = result_path + '/result/GetExp/tissueCut/segmentation/TissueCut.log'
		self.final_result_json = outfile_path + '/new_final_result.json'
		#self.saturation_file = result_path + "/result/GetExp/sequence_saturation.txt"
		self.result_dict = defaultdict(dict)
		self.result_dict['version'] = "version_v2"
	
		os.makedirs(self.outfile_path, exist_ok=True)
		

	def run(self):
		if os.path.exists(self.filter_result_path):
			self.filter_result()
		if os.path.exists(self.alignment_result_path):
			self.alignment_result()
#		if os.path.exists(self.exp_result_path):
#			self.basic()
		if os.path.exists(self.tissueCut):
			self.tissueCutStat()
		if os.path.exists(self.cell_cluster):
			self.copy2()
		if os.path.exists(self.figure_path):
			self.copy()
		#if os.path.exists(self.saturation_file):
		#	self.saturationPlot()
		with open(self.final_result_json,'w') as f:
			t_result_dict = json.loads(json.dumps(self.result_dict).replace('Umi','MID').replace('umi','MID').replace('barcode','CID').replace('Barcode','CID'))
			j = json.dump(t_result_dict, f, indent=4)
#			j = json.dump(self.result_dict, f, indent=4)

	def number_transfer(self, num):
		if isinstance(num,int) or isinstance(num,float):
			return(self.change(num))
		elif isinstance(num,str):
			if num.endswith('%'):
				return(num.replace(" ",""))
			else:
				return(self.change(num))
		elif isinstance(num,list):
			number_list = num
			if len(number_list) == 1:
				if number_list[0].endswith('%'):
					return(number_list[0])
				else:
					return(self.change(number_list[0]))
			#elif len(number_list) == 2:
			#	if number_list[1].endswith('%'):
			#		return(number_list[1])
			#	else:
			#		return(self.change(number_list[1]))
			elif len(number_list) == 2:
				number = number_list[0]
				percentage = number_list[1]
				if percentage.endswith('%'):
					final_percentage = percentage
				else:
					final_percentage = str(round(float(percentage)*100,2))+'%'
				number_changed = self.change(number)
				return(str(number_changed)+'('+str(final_percentage)+')')

	def change(self, number_raw):
		number = float(number_raw)
		if number > 1000000000000:
			return(str(round(number/1000000000000,2))+'T')
		elif number >1000000000:
			return(str(round(number/1000000000,2))+'G')
		elif number >1000000:
			return(str(round(number/1000000,2))+'M')
		elif number >1000:
			return(str(round(number/1000,2))+'K')
		else:
			return(str(round(number,2)))

	def copy(self):
		cp_bin_figures = 'cp ' + self.figure_path + '* ' + self.outfile_path
		os.system(cp_bin_figures)

	def copy2(self):
		cp_cell_cluster = 'cp ' + self.cell_cluster + ' ' + self.outfile_path
		os.system(cp_cell_cluster)

	def tissueCutStat(self):
		tissueCut_dict = defaultdict(list)
		line_num = 0
		with open(self.tissueCut, 'r') as file_tissue_cut:
			total_dict = defaultdict()
			name_bin_cut_list = []
			value_bin_cut_list = []
			for line in file_tissue_cut:
				line_num += 1
				if line_num == 1:
					if re.search('Cell', line.strip()):
						sub_title = '4.3.CellCut_Bin_stat'
						total_title = '4.1.CellCut_Total_Stat'
						final_title = '4.CellCut'
					else:
						sub_title = '4.2.TissueCut_Bin_stat'
						total_title = '4.1.TissueCut_Total_Stat'
						final_title = '4.TissueCut'
				elif line_num < 9:
					name = line.strip().split(':')[0]
					value = self.number_transfer(line.strip().split(':')[1])
					total_dict[name] = value
				else:
					if re.search('=',line):
						name = line.strip().split('=')[0]
						value = line.strip().split('=')[1]
						name_bin_cut_list.append(name)
						value_bin_cut_list.append(value)
					elif re.search(':',line):
						name = line.strip().split(':')[0]
						value = self.number_transfer(line.strip().split(':')[1])
						name_bin_cut_list.append(name)
						value_bin_cut_list.append(value)
		bin_line_num = 0
		if sub_title == '4.2.TissueCut_Bin_stat':
			for i in range(len(name_bin_cut_list)):
				bin_line_num += 1
				if bin_line_num % 7 == 1:
					temp_bin_cut_dict = {}
					temp_bin_cut_dict[name_bin_cut_list[i]] = value_bin_cut_list[i]
				elif bin_line_num % 7 == 0:
					temp_bin_cut_dict[name_bin_cut_list[i]] = value_bin_cut_list[i]
					tissueCut_dict[sub_title].append(temp_bin_cut_dict)
				else:
					temp_bin_cut_dict[name_bin_cut_list[i]] = value_bin_cut_list[i]
		else:
			temp_bin_cut_dict = {}
			for i in range(len(name_bin_cut_list)):
				bin_line_num += 1
				if bin_line_num < 11:
					temp_bin_cut_dict[name_bin_cut_list[i]] = value_bin_cut_list[i]
			tissueCut_dict[sub_title].append(temp_bin_cut_dict)

		tissueCut_dict[total_title].append(total_dict)				

		self.result_dict[final_title] = tissueCut_dict

	def basic(self):
		basic_dict = defaultdict(list)
		thread_list = []
		for fnames in os.listdir(self.bin_path):
			if fnames.endswith('gene.txt'):
				command1 = 'cp '+ self.bin_path + fnames + ' '+ self.outfile_path
				os.system(command1)
				command_dnb = self.python + ' ' + self.get_bin_data_py + ' -i ' + self.bin_path + fnames + ' -o ' + self.outfile_path + '/dnb_merge -t dnb -s 1000 '
				t1 = myThread(command_dnb)
				thread_list.append(t1)
				command_gene = self.python + ' ' + self.get_bin_data_py + ' -i ' + self.bin_path + fnames + ' -o ' + self.outfile_path + '/gene_merge -t gene -s 1000 '
				t2 = myThread(command_gene)
				thread_list.append(t2)

		for t in thread_list:
			t.start()
	
		for t in thread_list:
			t.join()
				
	def filter_result(self):
		adapter_filter_dict = defaultdict(list)
		list_sample_list = []
		for fnames in os.listdir(self.filter_result_path):
			if fnames.endswith('_barcodeMap.stat'):
				list_sample_list.append(fnames)
		for list_index in range(len(list_sample_list)):
			stat_id = list_sample_list[list_index]
			real_sample_id = stat_id.split('_barcodeMap.stat')[0]
			stat_id = os.path.join(self.filter_result_path, stat_id)
			name_key_list = ['getBarcodePositionMap_uniqBarcodeTypes','total_reads','mapped_reads', 'reads_with_adapter', 'fixed_sequence_contianing_reads','barcode_misOverlap_reads','barcode_withN_reads','Q10_bases_in_barcode','Q20_bases_in_barcode','Q30_bases_in_barcode','Q10_bases_in_umi','Q20_bases_in_umi','Q30_bases_in_umi', 'Q10_bases_in_seq', 'Q20_bases_in_seq', 'Q30_bases_in_seq', 'umi_filter_reads','umi_with_N_reads','umi_with_polyA_reads','umi_with_low_quality_base_reads','reads_with_adapter', 'reads_with_dnb']
			q30_value_list = []
			total_base_list = []
			if os.path.exists(stat_id):
				temp_dict = {}
				filter_dict = {}
				map_dict = {}
				with open(stat_id,'r') as file_stat_id:
					for line_file_stat_id in file_stat_id:
						for name_key in name_key_list:
							if re.search(name_key,line_file_stat_id):
								map_dict['Sample_id'] = real_sample_id
								name = line_file_stat_id.strip().split(':')[0]
								value = line_file_stat_id.strip().split(':')[1].strip().split('\t')
								#value = self.number_transfer(number)
								if (name == "reads_with_adapter"):
									filter_dict['reads_with_adapter'] = value[0]
								elif (name == "reads_with_dnb"):
									filter_dict['reads_with_dnb'] = value[0]
								else:
									temp_dict[name] = value[0]
									map_dict[name] = self.number_transfer(value)
				filter_dict['Sample_Name'] = real_sample_id
				filter_dict['Raw_Reads'] = self.number_transfer(int(temp_dict['mapped_reads']) - int(temp_dict['umi_filter_reads']))
				filter_dict['Clean_Reads'] = self.number_transfer(int(temp_dict['mapped_reads']) - int(temp_dict['umi_filter_reads']) - int(filter_dict['reads_with_adapter']) - int(filter_dict['reads_with_dnb']))
				filter_dict['Low_Quality_Reads'] = '0'
				filter_dict['Too_Many_N_Reads'] = '0'
				filter_dict['Too_Long_Reads'] = '0'
				filter_dict['Remarks'] = '0'
				filter_dict['reads_with_adapter'] = self.number_transfer(filter_dict['reads_with_adapter'])
				filter_dict['reads_with_dnb'] = self.number_transfer(filter_dict['reads_with_dnb'])

				adapter_filter_dict['1.1.Adapter_Filter'].append(map_dict)
				adapter_filter_dict['1.2.Filter_Stat'].append(filter_dict)
						
		self.result_dict['1.Filter_and_Map'] = adapter_filter_dict

	def alignment_result(self):
		alignment_dict = defaultdict(list)
		list_sample_id_list = []

		temp_dedup_dict = {}
		temp_annotation_dict = {}

		for fnames in os.listdir(self.alignment_result_path):
			if fnames.endswith('.Log.final.out'):
				sample_id = fnames.split('.Log.final.out')[0]
				star_log = self.alignment_result_path + fnames
				line_num = 0
				temp_input_read_dict = {}
				temp_uniq_dict = {}
				temp_multi_dict = {}
				temp_unmapping_dict = {}
				temp_chimeric_dict ={}
				with open(star_log,'r') as file_star_log:
					for line_file_star_log in file_star_log:
						temp_input_read_dict['Sample_Id'] = sample_id
						temp_uniq_dict['Sample_Id'] = sample_id
						temp_multi_dict['Sample_Id'] = sample_id
						temp_unmapping_dict['Sample_Id'] = sample_id
						temp_chimeric_dict['Sample_Id'] = sample_id
						line_num = line_num + 1
						if line_num<5:
							continue
						elif line_num == 6:
							temp_input_read_dict['Number_Of_Input_Reads'] = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 7:
							temp_input_read_dict['Average_Input_Read_Length'] = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 9:
							temp_uniq_dict['Mapped_Reads_Number'] = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 10:
							temp_uniq_dict['Mapped_Reads(%)'] = line_file_star_log.strip().split('|')[1].replace("\t","")
						elif line_num == 11:
							temp_uniq_dict['Average_Mapped_Length'] = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 24:
							multiRead = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 25:
							multiRead_percentage = line_file_star_log.strip().split('|')[1].replace("\t","")
							temp_multi_dict['Multiple_Loci'] = multiRead + '(' + multiRead_percentage + ')'
						elif line_num == 26:
							manyRead = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 27:
							manyRead_percentage = line_file_star_log.strip().split('|')[1].replace("\t","")
							temp_multi_dict['Many_Loci'] = manyRead + '(' + manyRead_percentage + ')'
						elif line_num == 29:
							mismatch_read = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 30:
							mismatch_read_percentage = line_file_star_log.strip().split('|')[1].replace("\t","")
							temp_unmapping_dict['Too_Many_Mismatches'] = mismatch_read + '(' + mismatch_read_percentage + ')'
						elif line_num == 31:
							too_short_read = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 32:
							too_short_read_percentage = line_file_star_log.strip().split('|')[1].replace("\t","")
							temp_unmapping_dict['Too_Short'] = too_short_read + '(' + too_short_read_percentage + ')'
						elif line_num == 33:
							other_read = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 34:
							other_read_percentage = line_file_star_log.strip().split('|')[1].replace("\t","")
							temp_unmapping_dict['Other'] = other_read + '(' + other_read_percentage + ')'
						elif line_num == 36:
							chimeric_read = self.number_transfer(line_file_star_log.strip().split('|')[1].replace("\t",""))
						elif line_num == 37:
							chimeric_read_percentage = line_file_star_log.strip().split('|')[1].replace("\t","")
							temp_chimeric_dict['Number_Of_Chimeric_Reads'] = chimeric_read + '(' + chimeric_read_percentage + ')'
				alignment_dict['2.1.Input_Read'].append(temp_input_read_dict)
				alignment_dict['2.2.Uniquely_Mapped_Read'].append(temp_uniq_dict)
				alignment_dict['2.3.Multi_Mapping_Read'].append(temp_multi_dict)
				alignment_dict['2.4.Unmapping_Read'].append(temp_unmapping_dict)
				alignment_dict['2.5.Chimeric_Read'].append(temp_chimeric_dict)
			elif fnames.endswith('summary.stat'):
				dedup_file = self.alignment_result_path + fnames
				line_num = 0
				with open(dedup_file, 'r') as file_dedup_file:
					for line in file_dedup_file:
						line_num = line_num + 1
						if line_num == 2:
							name_dup_list = line.strip().split('\t')
						elif line_num == 3:
							value_dup_list = line.strip().split('\t')
							for i in range(len(name_dup_list)):
								if name_dup_list[i] in temp_dedup_dict.keys():
									temp_dedup_dict[name_dup_list[i]].append(float(value_dup_list[i]))
								else:
									temp_dedup_dict[name_dup_list[i]] = [float(value_dup_list[i])]
						elif line_num == 5:
							name_annotation_list = line.strip().split('\t')
						elif line_num == 6:
							value_annotation_list = line.strip().split('\t')
							for i in range(len(name_annotation_list)):
								if name_annotation_list[i] in temp_annotation_dict.keys():
									temp_annotation_dict[name_annotation_list[i]] = temp_annotation_dict[name_annotation_list[i]] + float(value_annotation_list[i])
								else:
									temp_annotation_dict[name_annotation_list[i]] = float(value_annotation_list[i])
		for k,v in temp_dedup_dict.items():
			if k != 'Sample_Id':
				if k == 'FAIL_FILTER_RATE' or k == 'FAIL_ANNOTATE_RATE' or k == 'DUPLICATION_RATE':
					temp_dedup_dict[k] = self.number_transfer(mean(v))
				else:
					temp_dedup_dict[k] = self.number_transfer(sum(v))

		for k,v in temp_annotation_dict.items():
			if k != 'Sample_Id':
				temp_annotation_dict[k] = self.number_transfer(v)

		alignment_dict['2.6.Filter_And_Deduplication'].append(temp_dedup_dict)
		alignment_dict['2.7.Annotation'].append(temp_annotation_dict)

		self.result_dict['2.Alignment'] = alignment_dict

	def saturationPlot(self):
		import matplotlib.pyplot as plt
		import pandas as pd
		os.makedirs(self.outfile_path, exist_ok=True)
		sadf = pd.read_csv(self.saturation_file, sep=" ")
		fig=plt.figure(figsize=(10,8),dpi=100)
		plt.clf()
		ax = plt.subplot(1, 2, 1)
		ax.plot(sadf['bin_x'], sadf['bar_y1'])
		ax.set_xlabel("Total reads number of sampling")
		ax.set_ylabel("Sequencing Saturation")
		ax.grid()
		ax = plt.subplot(1, 2, 2)
		ax.plot(sadf['bin_x'], sadf['bin_y2'])
		ax.set_xlabel("Total reads number of sampling")
		ax.set_ylabel("Median Genes per bin")
		ax.grid()
		figFilePath = os.path.join(self.outfile_path, "plot_150x150_saturation.png")
		plt.tight_layout()
		plt.savefig(figFilePath, format="png", bbox_inches='tight')

class myThread(threading.Thread):
	def __init__(self, command):
		threading.Thread.__init__(self)
		self.cmd = command

	def run(self):
#		print("Starting " + self.cmd)
		os.system(self.cmd)
#		print("Exiting " + self.cmd)


if __name__=="__main__":		
	main()
