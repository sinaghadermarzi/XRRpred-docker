import os
for i in range (2,1010):
	in_filename = './seqs/' + str(i)+'.seq'
	out_file_name  = './res_short/' + str(i)+'.res'
#	print("Running ProfBVAL on sequence" + in_filename)
	command = "./iupred "+in_filename+" short > " + out_file_name
	print(command)
	os.system(command)
