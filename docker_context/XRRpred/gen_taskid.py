from datetime import datetime
from random import choice
from string import digits

def gen_taskid(prog_name):
    random_code= (''.join(choice(digits) for i in range(8)))
    time_string=  datetime.now().strftime("%Y%m%d%H%M%S")
    taskid = time_string+"_"+prog_name+"_"+ random_code
    return taskid
