import subprocess
import os
iupred_cmd = "/home/sghadermarzi/IUPRED/iupred"



def pred_iupred(seq):
    # using module subprocess to pass arguments to and from an
    # external program:
    # put in the corresponding strings for
    # "mycmd" --> external program name
    # "myarg" --> arg for external program
    # several args can be in the list ["mycmd", "myarg1", "myarg2"]
    # might need to change bufsize
    #create a sequence file in the tmp dir in the current path
    if not os.path.exists('./tmp'):
        os.mkdir('./tmp')
    with open("./tmp/tmp_seq.fasta","w",newline = "") as seqfile:
        seqfile.writelines(">tmp\n"+seq)
    p = subprocess.Popen([iupred_cmd+" ./tmp/tmp_seq.fasta"+" long"], bufsize=2048, shell=True,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
    # allow external program to work
    p.wait()
    # read the result to a string
    result_str = str(p.stdout.read())
    res_list = list()
    lines  = result_str.split("\n")
    for line in lines[9:]:
        line = line.rstrip('\n')
        line_cols = line.split()
        value = (float(line_cols[2]))
        if value<0:
            print ('ERROR!')
        value = round(value,2)
        res_list.append(value)
    return result_str


sequence = "MAAEPVEDNCINFVAMKFIDNTLYFIAEDDENLESDYFGKLESKLSVIRNLNDQVLFIDQGNRPLFEDMTDSDCRDNAPRTIFIISMYKDSQPRGMAVTISVKCEKISTLSCENKIISFKEMNPPDNIKDTKSDIIFFQRSVPGHDNKMQFESSSYEGYFLACEKERDLFKLILKKEDEL"

print("IUPRED OUT:\n",pred_iupred(sequence))


