import subprocess
import os
import shutil
import sys
from datetime import datetime
from random import choice
from string import digits
from gen_taskid import gen_taskid

#folder paths without slash
DEBUG = True
this_script_directory = os.path.dirname(os.path.realpath(__file__))
iupred_dir = os.path.abspath(this_script_directory+"/iupred")
asaquick_dir = os.path.abspath(this_script_directory+"/asaquick/bin")

if not os.path.exists(this_script_directory+'/tmp'):
	os.mkdir(this_script_directory+'/tmp')
temp_dir = os.path.abspath(this_script_directory+'/tmp')



def iupred_long(seq):
	#save the current working directory
	cwd  = os.getcwd()
	taskid = gen_taskid("iupred_long")
	#go to iupred directory
	os.chdir(iupred_dir)
	#create a sequence file in the temp_dir
	with open(temp_dir+"/"+taskid+".fasta","w") as seqfile:
		seqfile.writelines(">tmp\n"+seq)
	p = os.system("./iupred "+temp_dir+"/"+taskid+".fasta long > "+temp_dir+"/"+taskid+".res")
	# allow external program to work
	# read the result to a string
	with open(temp_dir+"/"+taskid+".res") as res_file:
		result_str = res_file.read()
		res_list = list()
		lines  = result_str.split("\n")
		for line in lines[9:]:
			line = line.rstrip('\n')
			line_cols = line.split()
			if len(line)>2:
				value = (float(line_cols[2]))
				res_list.append(value)
	#remove temp fasta file
	if not DEBUG:
		os.remove(temp_dir+"/"+taskid+".fasta")
		os.remove(temp_dir+"/"+taskid+".res")
	#go back to the previous working directory
	os.chdir(cwd)
	return res_list
	


def iupred_short(seq):
	taskid = gen_taskid("iupred_short")
	#save the current working directory
	cwd  = os.getcwd()
	#go to iupred directory
	os.chdir(iupred_dir)
	#create a sequence file in the temp_dir
	with open(temp_dir+"/"+taskid+".fasta","w") as seqfile:
		seqfile.writelines(">tmp\n"+seq)
	p = os.system("./iupred "+temp_dir+"/"+taskid+".fasta short > "+temp_dir+"/"+taskid+".res")
	# allow external program to work
	# read the result to a string
	with open(temp_dir+"/"+taskid+".res") as res_file:
		result_str = res_file.read()
		res_list = list()
		lines  = result_str.split("\n")
		for line in lines[9:]:
			line = line.rstrip('\n')
			line_cols = line.split()
			if len(line)>2:
				value = (float(line_cols[2]))
				res_list.append(value)
	#remove temp fasta file
	if not DEBUG:
		os.remove(temp_dir+"/"+taskid+".fasta")
		os.remove(temp_dir+"/"+taskid+".res")
	#go back to the previous working directory
	os.chdir(cwd)
	return res_list



def asaquick(seq):
	#this is a wrapper function for the asaquick tool for prediction of solvent accessible area
	#it takes the input sequence directly as a string (only containing aminoacid sequence) and outputs the a list of numbers containing residue level predictions
	#the numbers outputed from the program are transformed to be between 0 and 1: 0 means fully buried and 1 means fully accessible

	#save the current working directory
	cwd  = os.getcwd()
	#add asaquick to PATH
	sys.path
	sys.path.append(asaquick_dir)
	# initialize temp folder: where temporary files are stored
	# this_script_directory = os.path.dirname(os.path.realpath(__file__))
	# if not os.path.exists(this_script_directory + '/tmp'):
	#     os.mkdir(this_script_directory + '/tmp')
	# temp_dir = this_script_directory + '/tmp'
	taskid = gen_taskid("asaquick")
	#create a sequence file in the temp_dir
	os.chdir(asaquick_dir)
	with open(temp_dir+"/"+taskid+".fasta","w") as seqfile:
		seqfile.writelines(">tmp\n"+seq)
	p = os.system("./ASAquick "+temp_dir+"/"+taskid+".fasta > /dev/null 2>&1")
	#read the results from ASAquick
	fileaddress = "./asaq."+taskid+".fasta/rasaq.pred"
	f = open(fileaddress)
	res_list = list()
	line = f.readline()
	while line:
		line = line.rstrip('\n')
		line_cols = line.split(" ")
		# value = float(line_cols[2])
		value = (float(line_cols[2])/2)+0.5
		res_list.append(value)
		# res_list.append(line_cols[2])
		# use realine() to read next line
		line = f.readline()
	res_list = res_list[:-1]
	#remove temp fasta file
	if not DEBUG:
		os.remove(temp_dir+"/"+taskid+".fasta")
		#remove the results from ASAquick
		shutil.rmtree("./asaq."+taskid+".fasta")
	#go back to the previous working directory
	os.chdir(cwd)
	return res_list




if (__name__ == "__main__"):
	#run this script if you want to test these functions
	test_seq = "MEAENAGSYSLQQAQAFYTFPFQQLMAEAPNMAVVNEQQMPEEVPAPAPAQEPVQEAPKGRKRKPRTTEPKQPVEPKKPVESKKSGKSAKSKEKQEKITDTFKVKRKVDRFNGVSEAELLTKTLPDILTFNLDIVIIGINPGLMAAYKGHHYPGPGNHFWKCLFMSGLSEVQLNHMDDHTLPGKYGIGFTNMVERTTPGSKDLSSKEFREGGRILVQKLQKYQPRIAVFNGKCIYEIFSKEVFGVKVKNLEFGLQPHKIPDTETLCYVMPSSSARCAQFPRAQDKVHYYIKLKDLRDQLKGIERNMDVQEVQYTFDLQLAQEDAKKMAVKEEKYDPGYEAAYGGAYGENPCSSEPCGFSSNGLIESVELRGESAFSGIPNGQWMTQSFTDQIPSFSNHCGTQEQEEESHA"
	print("iupred long output for example :"+test_seq)
	print("is : "+ str(iupred_long(test_seq)))
	print("iupred short output for example :"+test_seq)
	print("is : "+ str(iupred_short(test_seq)))
	print("asaquick output for example :"+test_seq)
	print("is : "+ str(asaquick(test_seq)))