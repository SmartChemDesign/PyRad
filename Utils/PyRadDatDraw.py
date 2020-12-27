import sys
import numpy as np
from tqdm import tqdm
from re import findall
from sys import argv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


ts = None
data = dict()


def read_dat(file):
    """
    Read .dat file
    :param file: path&filename
    :return: None
    """
    global data
    global ts
    
    # Deploy iteration reader routine
    bDeployIterRead = False
    
    # check the count of lines in the file
    num_lines = sum(1 for _ in open(file, 'rb'))
    
    # process the file!
    with open(file, 'r') as f:
        for line in tqdm(f, total=num_lines, mininterval=0.5):
            line = line.strip()
            
            # Read time step
            if (bDeployIterRead == False) and (line.startswith("used time step: ")):
                ts = float(findall(r"\d+\.\d+", line.split(':')[1])[0])
                continue
            
            # Read iterations
            if line.startswith("iteration "):
                bDeployIterRead = True
                continue
                
            if bDeployIterRead:
                k, v = line.split(':')
                v = float(v)
                if k not in data:
                    data[k] = list()
                data[k].append(v)
                continue

def draw_plots():
    """
    Draws a plot and saves an image file
    :return: None
    """
    matplotlib.rcParams['savefig.dpi'] = 600
    matplotlib.rcParams["figure.dpi"] = 100
    
    # plot configuration
    fig, ax = plt.subplots()
    
    ax.set_xlabel('time, fs')
    ax.set_ylabel('species count')

    x_axis = np.asarray(range(len(data["h2o"]))) * ts

    if sum(data["h2"]) > 0:
        ax_h2, = ax.plot(x_axis, data["h2"], marker='', color='y',
                     linewidth=1.5, label="$H_2$")
    if sum(data["h2o2"]) > 0:                     
        ax_h202, = ax.plot(x_axis, data["h2o2"], "-", color='m',
                       linewidth=2, label="$H_2O_2$")
    if sum(data["h_diss"]) > 0:                   
        ax_h_diss, = ax.plot(x_axis, data["h_diss"], ":", color='g',
                         linewidth=2.5, label="$H^*$")
                         
    if sum(data["o_diss"]) > 0:                     
        ax_o_diss, = ax.plot(x_axis, data["o_diss"], "--", color='r',
                         linewidth=1, label="$O^*$")
                         
    if sum(data["oh_diss"]) > 0:                     
        ax_oh_diss, = ax.plot(x_axis, data["oh_diss"], ":", color='b',
                          linewidth=1.5, label="$OH^*$")
                          
    if sum(data["h3o"]) > 0:                      
        ax_h3o, = ax.plot(x_axis, data["h3o"], marker='', color='c',
                      linewidth=1.5, label="$H_3O^*$")
                      
    if sum(data["ho2"]) > 0:                  
        ax_ho2, = ax.plot(x_axis, data["ho2"], marker='', color='k',
                      linewidth=1, label="$HO_2^*$")

    # use 2nd axis for water count
    ax_water = ax.twinx()
    ax_h2o, = ax_water.plot(x_axis, data["h2o"], ":", marker='', color='r',
                            linewidth=1, label="$H_2O$")
    ax_water.set_ylabel('$H_2O$ count', color="r")

    ax.legend(bbox_to_anchor=(1.13, 1))
    ax_water.legend(loc=2)

    # draw plot
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_water.yaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.show()
                        
if len(sys.argv) != 2:
    print("Usage: put .sys file as script arguments")
    exit()

read_dat(sys.argv[1])
draw_plots()
