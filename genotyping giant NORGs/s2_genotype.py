#!/usr/bin/python
#sv_geno

import sys,os
from collections import OrderedDict

script,sv_type,sv_refReads_file,sv_nonrefReads_file,sv_geno_out = sys.argv

def rd_num(sv_Reads_file):
	sv_Read_dict = OrderedDict()
	sv_name_all = []
	with open(sv_Reads_file) as inf:
		allsp = []
		for i in inf:
			lsp = i.split()
			if i.startswith("svName"):
				for sp in lsp[1:]:
					sv_Read_dict[sp] = OrderedDict()
					allsp.append(sp)
			else:
				sv_name = lsp[0]
				sv_name_all.append(sv_name)
				sp_num = len(lsp[1:])
				for index  in range(sp_num):
					rd_num = lsp[1:][index]
					total_rdNum = int(rd_num.split("#")[0]) + int(rd_num.split("#")[1])
					sv_Read_dict[allsp[index]] [sv_name] = total_rdNum
	return sv_Read_dict,sv_name_all



#0 indicating "no exists"; 1 indicating "hetero genotype";2 indicating "exits";un indicating "unknwon"
if __name__  == "__main__":
	geno_sv = OrderedDict()
	sv_refRead_dict_svList = rd_num(sv_refReads_file)
	sv_refRead_dict = sv_refRead_dict_svList[0]
	sv_name_lst_ref = sv_refRead_dict_svList[1]
	sv_nonrefRead_dict = rd_num(sv_nonrefReads_file)[0]
	for sp in sv_refRead_dict:
		geno_sv[sp] = OrderedDict()
		for sv in sv_refRead_dict[sp]:
			ref_num = sv_refRead_dict[sp][sv]
			nonref_num = sv_nonrefRead_dict[sp][sv]
			if ref_num > 0 and nonref_num == 0:
					geno_sv[sp][sv] = '0'
			elif ref_num == 0 and nonref_num > 0:
					geno_sv[sp][sv] = '2'
			elif ref_num > 0 and nonref_num > 0:
				if (ref_num/float(nonref_num) >= 3 and nonref_num <= 3):
					geno_sv[sp][sv] = '0'
				elif (ref_num/float(nonref_num) <= float(1)/3 and ref_num <= 3):
					geno_sv[sp][sv] = '2'
				elif ref_num <= 2 and nonref_num <= 2:
					geno_sv[sp][sv] = 'un'
				else:
					geno_sv[sp][sv] = '1'
			else:
				geno_sv[sp][sv] = 'un'
	with open(sv_geno_out,'w') as outf:
		print >>outf,'\t'.join(['sp_name'] + sv_name_lst_ref)
		for sp in geno_sv:
			geno_sp_sv = []
			for sv in geno_sv[sp]:
				geno_sp_sv.append(geno_sv[sp][sv])
			print >>outf,'\t'.join([sp]+geno_sp_sv)
			
