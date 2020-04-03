#!/usr/bin/python
import pysam
import os 
import sys
from glob import glob
from collections import OrderedDict

#fetch(,,):Note `start` and `stop` denote 0-based
script,svfile,svType,bamType,bamfilePath = sys.argv

def handleSam(samfile,chr,start,end,correspond_svSize):
	readCovNum = [0,0]
	readTotal = [0,0]
	readPairNum = [0,0]
	if start < end:
		for read in samfile.fetch(chr,start-1,end):
			read_s = read.reference_start+1
			read_e = read.reference_end
			rdn = read.query_name+'_1' if read.is_read1 else read.query_name+'_2'
			rdn1 = read.query_name+'_2' if read.is_read1 else read.query_name+'_1' 
			read_len = read.infer_read_length()
			n_read_s = read.next_reference_start+1
			# read_name = []
			read_name_left = []
			read_name_right = []
			#stat pair read numbers
			if read_s < start-3 and read_e > start+3:
			#if read_e - start >= 5 and read_s < start-2:
				# print read.query_name
				readCovNum[0] = readCovNum[0]+1
				if rdn not in read_name_left:
					read_name_left.append(rdn)
					readTotal[0] = readTotal[0] +1
			#stat right breakpoint cover reads #v4  end+2
			elif read_s < end-3 and read_e > end+3:
			#elif end - read_s >= 5 and read_e > end+2:
				readCovNum[1] = readCovNum[1]+1
				if rdn not in read_name_right:
					read_name_right.append(rdn)
					readTotal[1] = readTotal[1] +1
			elif read_s  >= start -3 and read_e <= end-3:   #version 3 
				if read.is_proper_pair == False:
					pass
				else:
					if n_read_s < read_s: 
						if (read_s-n_read_s+1+read_len) <= 600:  # version5:insertsize size should be smaller than 600 bp
							if n_read_s  < start-3:
								readPairNum[0] = readPairNum[0] +1
								if rdn not in read_name_left:
									read_name_left.append(rdn)
									readTotal[0] = readTotal[0] +1
							else:
								pass
					else:
						if abs(n_read_s-read_s+1+read_len) <= 600:  #insertsize size should be smaller than 600 bp
							if n_read_s > end +3:
								readPairNum[1] = readPairNum[1] + 1
								if rdn not in read_name_right:
									read_name_right.append(rdn)
									readTotal[1] = readTotal[1] +1
			elif read_e < start-3 and n_read_s > end+3:
				if correspond_svSize > 600 and (n_read_s-read_s+1+read_len) <= 600:
					readPairNum[0] = readPairNum[0] +1
					readPairNum[1] = readPairNum[1] + 1
					if rdn not in read_name_left and rdn not in  read_name_right:
						readTotal[0] = readTotal[0] +1
						readTotal[1] = readTotal[1] +1
						read_name_left.append(rdn)
						read_name_right.append(rdn)
			elif n_read_s < start -3 and read_s > end+3:
				if correspond_svSize > 600 and (read_s-n_read_s+1+read_len) <= 600:
					if rdn not in read_name_left and rdn not in  read_name_right:
						readPairNum[0] = readPairNum[0] +1
						readPairNum[1] = readPairNum[1] + 1
						readTotal[0] = readTotal[0] +1
						readTotal[1] = readTotal[1] +1
						read_name_left.append(rdn)
						read_name_right.append(rdn)
			else:
				pass
	else:
		#stat cover read numbers and pair reads
		for read in samfile.fetch(chr,end-1,start):
			read_s = read.reference_start+1
			read_e = read.reference_end   
			read_len = read.infer_read_length()
			n_read_s = read.next_reference_start+1
			rdn = read.query_name+'_1' if read.is_read1 else read.query_name+'_2'
			rdn1 = read.query_name+'_2' if read.is_read1 else read.query_name+'_1' 
			read_name_left = []
			read_name_right = []
			if read_s < end-3 and  read_e > start +3:
				readCovNum[0] += 1
				readCovNum[1] += 1
				if rdn not in read_name_left and rdn not in  read_name_right:
					readTotal[0] = readTotal[0] +1
					readTotal[1] = readTotal[1] +1
					read_name_left.append(rdn)
					read_name_right.append(rdn)
			elif  read_s < end-3 and n_read_s > start:
				if  correspond_svSize > 600 and (n_read_s-read_s+1+read_len) <= 600:
					readPairNum[0] = readPairNum[0] +1
					readPairNum[1] = readPairNum[1] + 1
					if rdn not in read_name_left and rdn not in  read_name_right:
						readTotal[0] = readTotal[0] +1
						readTotal[1] = readTotal[1] +1
						read_name_left.append(rdn)
						read_name_right.append(rdn)
			elif n_read_s < end-3 and read_e > start+3:
				if  correspond_svSize > 600 and (read_s-n_read_s+1+read_len) <= 600:
					readPairNum[0] = readPairNum[0] +1
					readPairNum[1] = readPairNum[1] + 1
					if rdn not in read_name_left and rdn not in  read_name_right:
						readTotal[0] = readTotal[0] +1
						readTotal[1] = readTotal[1] +1
						read_name_left.append(rdn)
						read_name_right.append(rdn)
			else:
				pass
	# readSE = [readCovNum[0]+readPairNum[0],readCovNum[1]+readPairNum[1]]
	return readTotal,readPairNum


result_file = open('{}_{}_total_read.txt'.format(svType,bamType),'w')
result_file01 = open('{}_{}_pair_read.txt'.format(svType,bamType),'w')
# svType: del and ins
# bamType: ref and nonref
bamlist = glob("{}/*bam".format(bamfilePath))
sp_list = [os.path.split(i)[-1].split('_')[0] for i in bamlist]
print >>result_file,'\t'.join(['svName'] + sp_list)
print >>result_file01,'\t'.join(['svName'] + sp_list)

sv_sr = OrderedDict()
sv_sr_pair = OrderedDict()

for bam in bamlist:
	samfile = pysam.AlignmentFile(bam, "rb")
	if (svType == "del" and bamType == "ref") or (svType == "ins" and bamType == "ref"):
		with open(svfile) as inf:
			for line in inf:
				lsp = line.split()
				sv_c = lsp[0]
				sv_s = int(lsp[1])
				sv_e = int(lsp[2])
				correspond_s = int(lsp[7])
				correspond_e = int(lsp[8])
				correspond_svSize = correspond_e-correspond_s+1 if correspond_s  < correspond_e else abs(correspond_e-correspond_s+1)
				svName = lsp[9]
				# sv_rdNum = [] + svName
				if svName not in sv_sr and svName not in sv_sr_pair:
					sv_sr[svName] = svName
					sv_sr_pair[svName] = svName
					read_sp = handleSam(samfile,sv_c,sv_s,sv_e,correspond_svSize)
					sv_sr[svName] = sv_sr[svName]+'\t{}#{}'.format(read_sp[0][0],read_sp[0][1])
					sv_sr_pair[svName] = sv_sr_pair[svName]+'\t{}#{}'.format(read_sp[1][0],read_sp[1][1])
				else:
					read_sp = handleSam(samfile,sv_c,sv_s,sv_e,correspond_svSize)
					sv_sr[svName] = sv_sr[svName]+'\t{}#{}'.format(read_sp[0][0],read_sp[0][1])
					sv_sr_pair[svName] = sv_sr_pair[svName]+'\t{}#{}'.format(read_sp[1][0],read_sp[1][1])
	elif (svType == "del" and bamType == "nonref") or (svType == "ins" and bamType == "nonref"):
		with open(svfile) as inf:
			for line in inf:
				lsp = line.split()
				sv_c = lsp[6]
				sv_s = int(lsp[7])
				sv_e = int(lsp[8])
				svName = lsp[9]
				correspond_s = int(lsp[1])
				correspond_e = int(lsp[2])
				correspond_svSize = correspond_e-correspond_s+1 if correspond_s  < correspond_e else abs(correspond_e-correspond_s+1)
				if svName not in sv_sr and svName not in sv_sr_pair:
					sv_sr_pair[svName] = svName
					sv_sr[svName] = svName
					read_sp = handleSam(samfile,sv_c,sv_s,sv_e,correspond_svSize)
					sv_sr[svName] = sv_sr[svName]+'\t{}#{}'.format(read_sp[0][0],read_sp[0][1])
					sv_sr_pair[svName] = sv_sr_pair[svName]+'\t{}#{}'.format(read_sp[1][0],read_sp[1][1])
				else:
					read_sp = handleSam(samfile,sv_c,sv_s,sv_e,correspond_svSize)
					sv_sr[svName] = sv_sr[svName]+'\t{}#{}'.format(read_sp[0][0],read_sp[0][1])
					sv_sr_pair[svName] = sv_sr_pair[svName]+'\t{}#{}'.format(read_sp[1][0],read_sp[1][1])
	else:
		pass


for sv_n in sv_sr:
	print >>result_file,sv_sr[sv_n]

result_file.close()

for sv_n in sv_sr_pair:
	print >>result_file01,sv_sr_pair[sv_n]

result_file01.close()
