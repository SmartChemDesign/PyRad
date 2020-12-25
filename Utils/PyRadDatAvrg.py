from tqdm import tqdm
from re import findall
from sys import argv

# header saving dictionary
header = dict()
# saving header for the first time
first_header = True
# used simulation time step
ts = None
# list for all data
alldata = list()
# count of files
filecnt = len(argv)-1


def read_dat(file):
    """
    Read .dat file
    :param file: path&filename
    :return: list of dicts of all the species
    """
    global header
    global first_header
    global ts
    
    # Deploy dictionaries' similarity check
    bDeployDictCheck = False
    # Deploy iteration reader routine
    bDeployIterRead = False
    
    # all collected data
    file_list = list()
    # data of only one iter
    file_dict = dict()
    
    # check the count of lines in the file
    num_lines = sum(1 for _ in open(file, 'rb'))
    
    # process the file!
    with open(file, 'r') as f:
        for line in tqdm(f, total=num_lines, mininterval=0.5):
            line = line.strip()
            
            if line.startswith("iteration "):
                bDeployIterRead = True
                file_list.append(file_dict.copy())
                file_dict.clear()
                continue
                
            if bDeployIterRead:
                header_key, header_value = line.split(':')
                file_dict[header_key] = header_value
                continue
            
            # Start header-reading sequence
            if line.startswith("used distances in angstroms:"):
                bDeployDictCheck = True
                continue
            
            # Read time step
            if line.startswith("used time step: "):
                first_header = False
                bDeployDictCheck = False
                header_time = float(findall(r"\d+\.\d+", line.split(':')[1])[0])
                if ts is None:
                    ts = header_time
                elif ts != header_time:
                    # TODO: process files even if time steps are different
                    print("\nERROR: time steps don't match!")
                    exit()
                continue
            
            # Save distance string to header dictionary
            if bDeployDictCheck:
                header_key, header_value = line.split(':')
                if header_key in header:
                    # Check if file header is the same for all used files
                    if header[header_key] != header_value:
                        # Averaging of different kinds of simulation is useless
                        print("\nERROR: headers don't match!")
                        exit()           
                elif first_header:
                    # Save new header value
                    header[header_key] = header_value
                else:
                    # Again - usless task, quit
                    print("\nERROR: headers don't match!")
                    exit()                       
                continue
    file_list.append(file_dict.copy())
    file_dict.clear()
    return file_list.copy()


def find_avrg():
    """
    Finds average simulation
    :return: avrg simulation list
    """
    global alldata
    global filecnt
    
    avrgList = list()
    avrgDict = dict()
    
    # check all data have the same size   
    data_size = None
    for sim in alldata:
        size = len(sim)
        if data_size != size:
            if data_size is None:
                data_size = size
            else:
                print("ERROR: dataset lengths don't match!")
                exit()
   
    print("Summing...")
    # summing iteration over steps
    for i in tqdm(range(data_size)):
        avrgDict.clear()
        for sim in alldata:
            for k, v in sim[i].items():
                if k not in avrgDict:
                    avrgDict[k] = 0
                avrgDict[k] += int(v)
        avrgList.append(avrgDict.copy())

    print("Dividing...")
    # dividing every value by number of files
    for i in tqdm(range(data_size)):
        for k in avrgList[i].keys():
            avrgList[i][k] /= filecnt
            
    return avrgList.copy()


def save_dat(res):
    """
    Save data to the .dat file
    :param res: averaged simulation list
    :return: None
    """
    global header
    global ts
    
    # write output
    with open("average.dat", "w+") as outfile:
        outfile.write("used distances in angstroms:\n")
        for h in header:
            outfile.write(f"\t{h}:{header[h]}\n")
        outfile.write(f"used time step: {ts} fs\n")
        outfile.write("\n")
        
        for i in tqdm(range(len(res)-1)):
            outfile.write(f"iteration #{i}\n")
            for k, v in res[i+1].items():
                outfile.write(f"\t{k}: {v}\n")


# MAIN CYCLE

# show help message
if filecnt < 2:
    print("Arguments: two or more PyRad .sys files; then this script will find average .sys")
    exit()

# read all files
for arg in argv[1:]:
    print(f"Reading file {arg}")
    alldata.append(read_dat(arg))
    
# find average simulation
print("Calculating average simulation")
results = find_avrg()

# save results to .sys file
print("Saving data...")
save_dat(results)
print("OK")
