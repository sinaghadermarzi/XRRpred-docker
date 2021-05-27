import pandas
import os
import numpy
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib
matplotlib.use('Agg')

def findnn_idx(arr, x):
    r = len(arr)-1
    l = 0
    while l <= r:
        mid = l + (r - l) // 2;
        # Check if x is present at mid
        if arr[mid] == x:
            return mid
            # If x is greater, ignore left half
        elif arr[mid] < x:
            l = mid + 1
        # If x is smaller, ignore right half
        else:
            r = mid - 1
    if l == len(arr):
        return l - 1
    elif r == -1:
        return 0
    else:
        dist_l = arr[l] - x
        dist_r = x - arr[r]
        if dist_l > dist_r:
            return r
        else:
            return l

is_sorted = lambda a: numpy.all(a[:-1] <= a[1:])
def prepare_plot(save_path,pop_vals,point_loc,annot_txt,xlabel,ylabel):
    # rc('text', usetex=True)
    if not is_sorted(pop_vals):
        pop_vals = sorted(pop_vals)
    fig = plt.figure(figsize=(8,4))
    # ax = sns.distplot(pop_vals,bins = 20,kde_kws=dict(gridsize=1000,bw =0.1,c = "k"))
    nice_top = numpy.percentile(pop_vals,90)
    nice_bottom = numpy.percentile(pop_vals,10)
    nice_range= nice_top-nice_bottom


    ax = sns.kdeplot(pop_vals,gridsize=200,color = "k")
    ax.set_xlabel(xlabel,fontsize = 14)
    ax.set_ylabel(ylabel,fontsize = 14)
    x_lim_low, x_lim_high = ax.get_xlim()
    all_range = x_lim_high - x_lim_low
    points = ax.get_lines()[0].get_data()
    # piece = len(points[0])/10

    # x_high_lim = max(nice_top,point_loc*1.1)
    # x_low_lim = min(nice_bottom,point_loc*0.8)
    # ax.xlim(x_low_lim ,x_high_lim)
    # ax.xlim(0,x_high_lim)
    y_lim_low, y_lim_high = ax.get_ylim()
    # ax.set_ylim(y_lim_low, y_lim_high*1.3)

    divider_locs = [0] * 11
    percs = [0]*11
    divider_locs[0] = 0
    divider_locs[10] = len(points[0])-1
    percs[0] = min(pop_vals)
    percs[10] = max(pop_vals)


    for i in range(1,10):
        # percs[i] = numpy.percentile(points[0],i*10)
        percs[i] = numpy.percentile(pop_vals,i*10)
        divider_locs[i] = findnn_idx(points[0],percs[i])
        


    alpha = 0.8
    my_colors= [
        (0/255, 163/255, 0/255, alpha), # interval 1 (0-10 percentile)
        (0/255, 163/255, 0/255, alpha), # interval 2 (10-20 percentile)
        (145/255, 220/255, 145/255, alpha), # interval 3 (20-30 percentile)
        (145/255, 220/255, 145/255, alpha), # interval 4 (30-40 percentile)
        (200/255, 200/255, 200/255, alpha), # interval 5 (40-50 percentile)
        (200/255, 200/255, 200/255, alpha), # interval 6 (50-60 percentile)
        (220/255, 145/255, 145/255, alpha), # interval 1 (60-70 percentile)
        (220/255, 145/255, 145/255, alpha), # interval 1 (70-80 percentile)
        (255/255, 0/255, 0/255, alpha), # interval 1 (80-90 percentile)
        (255/255, 0/255, 0/255, alpha), # interval 1 (90-100 percentile)
    ]

    my_text_colors = [
        "palegreen", # interval 3 (20-30 percentile)
        "palegreen", # interval 4 (30-40 percentile)
        "green",# interval 1 (0-10 percentile)
        "green",# interval 2 (10-20 percentile)
        "olive", # interval 5 (40-50 percentile)
        "olive", # interval 6 (50-60 percentile)
        "darkred",# interval 1 (80-90 percentile)
        "darkred",# interval 1 (90-100 percentile)
        "darkred",
        "darkred"
    ]

    my_colors= [
        "green",# interval 1 (0-10 percentile)
        "green",# interval 2 (10-20 percentile)
        "palegreen", # interval 3 (20-30 percentile)
        "palegreen", # interval 4 (30-40 percentile)
        "palegoldenrod", # interval 5 (40-50 percentile)
        "palegoldenrod", # interval 6 (50-60 percentile)
        "lightcoral",# interval 1 (60-70 percentile)
        "lightcoral",# interval 1 (70-80 percentile)
        "tab:red",# interval 1 (80-90 percentile)
        "tab:red"# interval 1 (90-100 percentile)
    ]


    names =[
        "Very Good",
        "Very Good",
        "Good",
        "Good",        
        "Medium",
        "Medium",
        "Poor",
        "Poor",
        "Very Poor",
        "Very Poor"
    ]



    for i in range(len(my_colors)):
        c= my_colors[i]
        arr = list(colors.to_rgba(c))
        arr[3] = alpha
        my_colors[i] = tuple(arr)


    

    for i in range(0,5):
        idx1 = divider_locs[2*i]
        idx2= divider_locs[2*(i+1)]+1
        # mid = int((idx1+idx2)/2) 
        x = points[0][idx1:idx2]
        if i == 0:
            mid_x = (1.5*points[0][idx1]+8.5*points[0][idx2-1])/10
        elif i==4:
            mid_x = points[0][idx1]
        else:
            mid_x = (19*points[0][idx1]+points[0][idx2-1])/20
        # idxpoint = findnn_idx(x,point_loc)
        y = points[1][idx1:idx2]
        p = ax.fill_between(x,[0]*len(y),y,color=my_colors[2*i])
        t = ax.text(mid_x,y_lim_high/50,names[2*i],rotation="vertical",color = my_text_colors[2*i],fontsize=12)
    # ax.arrow(point_loc,y_lim_high*0.27,0,-y_lim_high*0.26,color  = (0,0,0,1),length_includes_head = True,head_length = y_lim_high/20,head_width=all_range/60)
    ax.annotate("", xy=(point_loc, 0), xytext=(point_loc, y_lim_high*0.26), arrowprops=dict(headlength=16,headwidth = 16,linewidth=1.5,facecolor=(0,0,0,1)))
    # , arrowprops={'arrowstyle': '->'}
    tb_x = point_loc-all_range*0.05
    tb_y = y_lim_high*0.45
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.70)
    ax.text(tb_x,tb_y,annot_txt, fontsize=10, verticalalignment='top', bbox=props)
    ax.tick_params(axis='both', which='major', labelsize=14)
    fig.savefig(save_path,transparent=True,bbox_inches='tight')
    return save_path

        
#         ax.fill(x_1,y_1,color = "None",hatch = "//",edgecolor ="cyan")
#         ax.set_edgecolor("cyan")
        
#         p = ax.fill_between(x,[0]*len(y),y,color=(0,0,0,0.99),hatch="/////",edgecolor = "cyan")
def q_label(perc):
	if perc<20:
		return "Very Good"
	elif perc>=20 and perc<40:
		return "Good"
	elif perc>=40 and perc<60:
		return "Medium"
	elif perc>=60 and perc<80:
		return "Poor"
	elif perc>=80:
		return "Very Poor"    

        
# for i in range(11):
#     line_x= points[0][line_locs[i-1]]
#     line_end_y = points[1][line_locs[i-1]]
#     ax.plot([line_x,line_x],[0,line_end_y],linewidth = 1,c = "k",linestyle='dashed')

def formatted(f):
	return format(f, '.3f').rstrip('0').rstrip('.')


def perc(pop,val):
	return int(100 * findnn_idx(pop,val)/len(pop))


def prepare_figure_tag(path_prefix, pid,pred_res,pred_rfree,current_profile_exit_code,current_profile_stop_chain,pop_res,pop_rfree):
	res_perc = perc(pop_res,pred_res)

	rfree_perc = perc(pop_rfree,pred_rfree)

	res_perc_str  =" (" + q_label(res_perc)+" Resolution; "+str(int(res_perc))+" Percentile)"
	res_annot_str =">"+pid+ "\nPredicted Resolution:\n"+formatted(pred_res)+" ("+q_label(res_perc)+"; "+str(int(res_perc))+r"%ile)"
	
	rfree_perc_str  = " (" + q_label(rfree_perc)+" R-Free; "+str(int(rfree_perc))+" Percentile)"
	rfree_annot_str =">"+ pid + "\nPredicted R-Free:\n" + formatted(pred_rfree) + " (" + q_label(rfree_perc) + "; " + str(int(rfree_perc)) + r"%ile)"


	# res_perc_str  = "" 
	# res_annot_str = ""
	# rfree_perc_str  = ""
	# rfree_annot_str = ""


    
	prepare_plot(path_prefix+"res_"+pid+".svg",pop_res,pred_res,res_annot_str,"Resolution (Angstrom)", "Probability Density")
	prepare_plot(path_prefix+"rfree_"+pid+".svg",pop_rfree,pred_rfree,rfree_annot_str,"R-Free", "Probability Density")

	




	tag_str= '<tr>\n'
	tag_str+= '<th colspan="8" class="seq" style="background-color:#FFA500; font-size: 18px;">>'+pid+'</th></tr>\n'
	if current_profile_exit_code=="success":
		tag_str+= '<tr><th colspan="8">\nPredicted Resolution:&nbsp'+formatted(pred_res)+'\n'
		tag_str+= res_perc_str+'<br/>\n'
		tag_str+= '<img style="display:block;" src="res_'+pid+'.svg" width="100%" ><br/>\n'
		tag_str+= 'Predicted R-Free:&nbsp'+formatted(pred_rfree)+'\n'
		tag_str+= rfree_perc_str+'<br/>\n'
		tag_str+= '<img style="display:block;" src="rfree_'+pid+'.svg" width="100%" ></th></tr>\n'
	else:
		error_text  = "Error in calculation of "+current_profile_exit_code+" for chain "+str(current_profile_stop_chain)
		tag_str+= '<tr><th colspan="8">'+error_text+'</th></tr>\n'
	return tag_str












def prepare_output(html_path, csv_path ,pid_list, res_list, rfree_list, profile_exit_code,profile_stop_chain, pop_res_sorted,pop_rfree_sorted):
    this_script_directory = os.path.dirname(os.path.realpath(__file__))
    os.chdir(this_script_directory)
    csv_path = os.path.abspath(csv_path)
    html_path = os.path.abspath(html_path)
    spl = csv_path.split("/")
    csv_address_suffix = "/".join(spl[3:])
    html_begining = r'''
    <!DOCTYPE html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta http-equiv="X-UA-Compatible" content="IE=Edge">
    <title>XRRPred Server - Results Page</title>
    <link rel="stylesheet" href="../../../css/bootstrap.css">
    <link rel="stylesheet" href="../../../assets/css/custom.min.css">
    <link rel="stylesheet" href="../../../servers/biomine.css">
    <link rel="stylesheet" href="../../stylesSCRIBER.css">
    </head>
    <body>
    <div class="container">
    <div class="row">
    <div class="col-lg-8 col-lg-offset-2">
    <h1>XRRPred results page</h1>
    <p>Results for <a target="_blank" href="http://biomine.cs.vcu.edu/servers/XRRPred/">XRRPred</a> webserver.</p>
    <p>Use this link to download the results as a CSV file: 
    <a target="_blank" href="http://biomine.cs.vcu.edu/'''+ csv_address_suffix +r'''">results.csv</a></p>
    </div>
    </div>
    <div class="row">
    <div class="col-lg-8 col-lg-offset-2">
    <div class="Predictions">
    <h2>Results</h2>
    <div class="table-responsive">
    <table class="table table-condensed">
    '''
    html_end  = r'''
    </table>
    </div>
    </div>
    </div>
    <div class="row">
    <div class="col-lg-8 col-lg-offset-2">
    <footer>
    <h2>Visit biomine lab web page</h2>
    <a target="_blank" href="http://biomine.cs.vcu.edu">http://biomine.cs.vcu.edu</a>
    </footer>
    </div>
    </div>
    </div>
    </body>
    </html>
    '''

    path_prefix = os.path.dirname(html_path)+"/"
    # with open("sdfsdfsdfsdf_path_prefix.txt","w") as cmdfile:
    #     cmdfile.writelines(path_prefix)
    with open(html_path ,"w") as htmlfile:
        htmlfile.writelines(html_begining)
        for i in range(len(pid_list)):
            htmlfile.writelines(prepare_figure_tag(path_prefix,pid_list[i],res_list[i],rfree_list[i],profile_exit_code[i],profile_stop_chain[i],pop_res_sorted,pop_rfree_sorted))
        htmlfile.writelines(html_end)

    df = pandas.DataFrame()
    df["Protein ID"]= pid_list
    df["Resolution Predicted by XRRPred (Angstrom)"]  = [formatted(x) for x in res_list]
    df["Percentile for Predicted Resolution"] =res_prec_list =  [perc(pop_res_sorted, x) for x in res_list]
    df["Qualitative Label for Predicted Resolution"] = [q_label(x) for x in res_prec_list ]
    df["R-Free Predicted by XRRPred"] = [formatted(x) for x in rfree_list]
    df["Percentile for Predicted R-Free"] = rfree_perc_list =[perc(pop_rfree_sorted, x) for x in rfree_list]
    df["Qualitative Label for Predicted R-Free"] = [q_label(x) for x in rfree_perc_list]

    for idx in range(len(pid_list)):
    	if profile_exit_code[idx] != "success":
    		df.iloc[idx,2:]=""
    		error_text  = "Error in calculation of "+profile_exit_code[idx]+" for chain "+str(profile_stop_chain[idx])



    df.to_csv(csv_path,index = False)


 


    # df= pandas.read_csv("training_set.csv")
    # res_vals = sorted(df["OUT_resolution"].values)
    # rfree_vals = sorted(df["OUT_rfree"].values)