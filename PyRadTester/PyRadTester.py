import os
import re
from operator import itemgetter, sub
from sys import argv

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from matplotlib.ticker import MaxNLocator
from matplotlib.widgets import Slider, Button
from scipy.spatial import distance


def get_system(file):
    """
    Loads whole system to the script
    No Element Tree anymore: sometimes we want to parse a broken file
    :param file: filename of the QB@LL .o file
    :return:
    1)ionic steps list  of atom list of dictionary (N, nom, coo)
    where coo is list of x,y,z coordinated
    2)step time in femtoseconds
    """
    # data for 1 ionic step
    ionic_step_data = list()
    # summary for all ionic steps
    ionic_step_summary = list()
    # summary for econst parameter
    econst_data = list()

    # compile regexp for separating atom symbol and its number
    pattern = re.compile("([a-zA-Z]+)([0-9]+)")

    # ionic step counter - parser's decision for ionic step incrementing
    stepnum = 1

    # flag to skip projectile
    skip = False

    # atom's name and number 
    num = 0
    name = ''

    # step time in atomic units
    au_time = 0.0

    num_lines = sum(1 for _ in open(file, 'rb'))

    with open(file, 'r') as f:
        for line in tqdm(f, total=num_lines, mininterval=0.5):
            line = line.rstrip("\n")
            # if line is set dt line - get au_time
            if line.startswith("<!-- [qball] set dt"):
                au_time = float(re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", line)[0])
                continue
                
            # if line is econst - add it to the list
            if line.startswith("  <econst>"):
                econst_data.append(float(re.findall(r"[-+]?\d*\.?\d+|[-+]?\d+", line)[0]))
                continue

            # if line is related to starting new ionic iteration
            if line.startswith("<iteration count="):
                iternum = int(line.split("\"")[1])
                # update stepnum and save the data
                if iternum != stepnum:
                    ionic_step_summary.append(list(ionic_step_data))
                    ionic_step_data.clear()
                    stepnum = iternum
                continue

            # if line is atom name - parse it            
            if line.startswith("  <atom name"):
                line = line.split("\"")[1]
                # except is rized if line is projectile info - skip it and its coos
                try:
                    line = pattern.match(line).groups()
                except AttributeError:
                    skip = True
                    continue
                # or save atom symbol and it's number
                num = int(line[1])
                name = line[0]
                continue

            # uf line is coordinates line and last element wasn't projectile - save it
            if line.startswith("    <position> "):
                # projectile check
                if skip:
                    skip = False
                    continue
                # save xyz
                line = line.split(" ")[5:8]
                # bohr -> angstroms
                coords = [float(coo) * 0.529177 for coo in line]
                # save dict to list
                ionic_step_data.append({"N": num, "nom": name, "coo": np.asarray(coords)})

        ionic_step_summary.append(list(ionic_step_data))

    # sort atoms by number, not by element symbol as QB@LL does
    for i in range(len(ionic_step_summary)):
        ionic_step_summary[i] = sorted(ionic_step_summary[i], key=itemgetter("N"))
        
    # save econst list to file
    np.savetxt(f"{name.split('.')[0] + '_econst.csv'}", np.asarray(econst_data), delimiter=',')
    
    return ionic_step_summary, au_time * 0.02418884254  # a.u. -> fs


def shape_h2o_mols_and_find_dist(system):
    """
    Modifies the same system list, but for oxygen there's additional list
    of h1 and h2 distances (H1 - O - H2)
    :param system: system to analyze
    :return: None
    """
    # loop by every system frame + corresponding number
    for frame in tqdm(system, mininterval=0.5):
        # compute distance between hydrogens and oxygen
        for i in range(1, len(frame), 3):
            h1 = distance.euclidean(frame[i]["coo"], frame[i - 1]["coo"])
            h2 = distance.euclidean(frame[i]["coo"], frame[i + 1]["coo"])
            frame[i]["dist"] = [h1, h2]
            

def count_h2o_mols_and_radicals(system, h_diss_dist, h2_creation_dist, h3o_creation_dist, ho2_creation_dist,
                                h2o2_creation_dist):
    """
    Returns only deviated molecules defined by criteria
    :param system: where to search
    :param h_diss_dist: distance when h counts as radical
    :param h2_creation_dist: distance when two h counts as h2 formation
    :param h3o_creation_dist: distance when water and H* creates H3O*
    :param ho2_creation_dist: distance when OH* and O* creates HO2*
    :param h2o2_creation_dist: distance between HO - OH when they create H2O2
    :return: list of different molecular species in the system
    """
    species_list = list()

    # loop by every frame in the system
    for frame in tqdm(system, mininterval=0.5):
        h2o = list()
        h_diss = list()
        oh_diss = list()
        o_diss = list()
        h2 = list()
        h3o = list()
        ho2 = list()
        h2o2 = list()

        # get info from oxygens (distances from O to H1 and H2)
        # 'primary' particles: H, OH
        for i in range(1, len(frame), 3):
            # H2O explosion!
            if (frame[i]["dist"][0] > h_diss_dist) and (frame[i]["dist"][1] > h_diss_dist):
                h_diss.append(frame[i - 1])
                h_diss.append(frame[i + 1])
                o_diss.append(frame[i])
                continue

            # dissociation OH + H
            if (frame[i]["dist"][0] > h_diss_dist) and not (frame[i]["dist"][1] > h_diss_dist):
                h_diss.append(frame[i - 1])
                oh_diss.append((frame[i], frame[i + 1]))
                continue

            # dissociation OH + H with another H
            if not (frame[i]["dist"][0] > h_diss_dist) and (frame[i]["dist"][1] > h_diss_dist):
                h_diss.append(frame[i + 1])
                oh_diss.append((frame[i], frame[i - 1]))
                continue

            # no dissociation
            h2o.append((frame[i - 1], frame[i], frame[i + 1]))

        # 'secondary' particles:

        # formation of H2 - iterate over all elements over all elements
        for item1 in h_diss:
            for item2 in h_diss:
                if item1 is not item2:
                    if distance.euclidean(item1["coo"], item2["coo"]) < h2_creation_dist:
                        h2.append((item1, item2))
                        # so item1 and item2 are not h_diss - remove them
                        h_diss.remove(item1)
                        h_diss.remove(item2)
                        break

        # formation of H3O - iterate over all water molecules and iterate over all h_diss
        for item1 in h2o:
            for item2 in h_diss:
                # item1[1] is oxygen atom
                if distance.euclidean(item1[1]["coo"], item2["coo"]) < h3o_creation_dist:
                    h3o.append((item1, item2))
                    h2o.remove(item1)
                    h_diss.remove(item2)
                    break

        # formation of HO2 as a product of OH_diss and O_diss
        for item1 in oh_diss:
            for item2 in o_diss:
                # item1[0] is oxygen atom
                if distance.euclidean(item1[0]["coo"], item2["coo"]) < ho2_creation_dist:
                    ho2.append((item1, item2))
                    # remove item1 and item 2 from old lists
                    oh_diss.remove(item1)
                    o_diss.remove(item2)
                    break

        # there's two places where we can find h2o2:
        # 1st - HO - OH
        for item1 in oh_diss:
            for item2 in oh_diss:
                if item1 is not item2:
                    # item1[0] and item2[0] are oxygen atoms
                    if distance.euclidean(item1[0]["coo"], item2[0]["coo"]) < h2o2_creation_dist:
                        h2o2.append((item1, item2, "HO-OH"))
                        # remove item1 and item 2 from old lists
                        oh_diss.remove(item1)
                        oh_diss.remove(item2)
                        break

        # 2nd - HO2 - H
        for item1 in ho2:
            for item2 in h_diss:
                # item1[1] is oxygen without hydrogen, H-O-(O) <- this one
                if distance.euclidean(item1[1]["coo"], item2["coo"]) < h_diss_dist:
                    h2o2.append((item1, item2, "HO2-H"))
                    # remove item1 and item 2 from old lists
                    ho2.remove(item1)
                    h_diss.remove(item2)
                    break
                    
        # re-creating of water molecules
        for item1 in h_diss:
            for item2 in oh_diss:
                # item2[0] is oxygen atom
                if distance.euclidean(item1["coo"], item2[0]["coo"]) < h_diss_dist:
                    h_diss.remove(item1)
                    oh_diss.remove(item2)
                    h2o.append((item1, item2[0], item2[1]))
                    break

        species_list.append({"h2o": h2o, "h2": h2, "h_diss": h_diss, "o_diss": o_diss, "oh_diss": oh_diss, "h3o": h3o,
                             "ho2": ho2, "h2o2": h2o2})
    return species_list


if __name__ == '__main__':
    # MAIN CYCLE========================================================================================================
    # TODO: rework matplotlib window. Now the best results only in 1920x1080 fullscreen

    # default distances - for sliders in matplotlib
    _h_diss_dist = 1.5
    _h2_creation_dist = 0.8
    _h3o_creation_dist = 2.9
    _ho2_creation_dist = 1.5
    _h2o2_creation_dist = 1.2

    # global names
    _name = ''
    species = None
    s = None

    def save_results(_):
        """
        Save button routine, save all species to .dat file
        :return: None
        """
        print("\tsaving dat file")
        global _name
        global species
        global dt
        # write output
        with open(_name.split(".")[0] + ".dat", "w+") as outfile:
            outfile.write(f"used distances in angstroms:\n\tH* dissociation: {_h_diss_dist}"
                          f"\n\tH2 creation: {_h2_creation_dist}\n\tH3O* creation: {_h3o_creation_dist}"
                          f"\n\tHO2* creation: {_ho2_creation_dist}\n\tH2O2 creation: {_h2o2_creation_dist}\n")
            outfile.write(f"used time step: {dt} fs\n\n")
            for i, item in enumerate(species):
                outfile.write(f"iteration #{i}\n")
                outfile.write(f"\th2o: {len(item['h2o'])}\n")
                outfile.write(f"\th_diss: {len(item['h_diss'])}\n")
                outfile.write(f"\to_diss: {len(item['o_diss'])}\n")
                outfile.write(f"\toh_diss: {len(item['oh_diss'])}\n")
                outfile.write(f"\th2: {len(item['h2'])}\n")
                outfile.write(f"\th3o: {len(item['h3o'])}\n")
                outfile.write(f"\tho2: {len(item['ho2'])}\n")
                outfile.write(f"\th2o2: {len(item['h2o2'])}\n")
        print("\tOK")


    def save_image(_):
        print("\tsaving jpg file")
        plt.savefig(f"{_name.split('.')[0]}.jpg", dpi=600,
                    bbox_inches=ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).expanded(1.8, 1.4))
                    

    def save_xyz(_):
        print("\tsaving xyz file")
        global _name
        global s
        global dt
        with open(_name.split(".")[0] + ".xyz", "w+") as xyzfile:
            for i, frame in enumerate(s):
                xyzfile.write(f"{len(frame)}\n")
                xyzfile.write(f"{dt*i} fs\n")
                for atom in frame:
                    xyzfile.write(f"{atom['nom']} {atom['coo'][0]} {atom['coo'][1]} {atom['coo'][2]}\n")
        print("\tOK")
        

    # check arguments line
    if len(argv) == 1:
        print("Arguments: .o qb@ll out file")
        exit()

    # plot configuration
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.45)
    ax.margins(x=0)

    axcolor = 'lightgoldenrodyellow'

    # create space for sliders
    ax_h_diss = plt.axes([0.25, 0.31, 0.65, 0.03], facecolor=axcolor)
    ax_h2 = plt.axes([0.25, 0.27, 0.65, 0.03], facecolor=axcolor)
    ax_h3o = plt.axes([0.25, 0.23, 0.65, 0.03], facecolor=axcolor)
    ax_ho2 = plt.axes([0.25, 0.19, 0.65, 0.03], facecolor=axcolor)
    ax_h2o2 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

    ax_save = plt.axes([0.8, 0.08, 0.12, 0.05])
    ax_saveimg = plt.axes([0.6, 0.08, 0.12, 0.05])
    ax_savexyz = plt.axes([0.4, 0.08, 0.12, 0.05])

    # create sliders
    S_h_diss = Slider(ax_h_diss, "H_diss", 0.1, 8.0, valinit=_h_diss_dist, valstep=0.1)
    S_h2 = Slider(ax_h2, "H2 creation", 0.1, 8.0, valinit=_h2_creation_dist, valstep=0.1)
    S_h3o = Slider(ax_h3o, "H3O creation", 0.1, 8.0, valinit=_h3o_creation_dist, valstep=0.1)
    S_ho2 = Slider(ax_ho2, "HO2 creation", 0.1, 8.0, valinit=_ho2_creation_dist, valstep=0.1)
    S_h2o2 = Slider(ax_h2o2, "H2O2 creation", 0.1, 8.0, valinit=_h2o2_creation_dist, valstep=0.1)
    
    B_save = Button(ax_save, "Save .dat")
    B_saveimg = Button(ax_saveimg, "Save .jpg")
    B_savexyz = Button(ax_savexyz, "Save .xyz")

    # create buttons
    B_save.on_clicked(save_results)
    B_saveimg.on_clicked(save_image)
    B_savexyz.on_clicked(save_xyz)


    def update(_):
        """
        Slider updater routine
        """
        global _h_diss_dist
        global _h2_creation_dist
        global _h3o_creation_dist
        global _ho2_creation_dist
        global _h2o2_creation_dist
        global species
        global s
        _h_diss_dist = S_h_diss.val
        _h2_creation_dist = S_h2.val
        _h3o_creation_dist = S_h3o.val
        _ho2_creation_dist = S_ho2.val
        _h2o2_creation_dist = S_h2o2.val
        print("\tfinding all molecular/atomic structures with new criteria...")
        species = count_h2o_mols_and_radicals(s, _h_diss_dist, _h2_creation_dist, _h3o_creation_dist,
                                              _ho2_creation_dist, _h2o2_creation_dist)

        ax_h2o.set_ydata([len(x["h2o"]) for x in species])
        ax_h2.set_ydata([len(x["h2"]) for x in species])
        ax_h202.set_ydata([len(x["h2o2"]) for x in species])
        ax_h_diss.set_ydata([len(x["h_diss"]) for x in species])
        ax_o_diss.set_ydata([len(x["o_diss"]) for x in species])
        ax_oh_diss.set_ydata([len(x["oh_diss"]) for x in species])
        ax_h3o.set_ydata([len(x["h3o"]) for x in species])
        ax_ho2.set_ydata([len(x["ho2"]) for x in species])

        ax.relim()
        ax_water.relim()
        ax.autoscale_view()
        ax_water.autoscale_view()

        fig.canvas.draw_idle()


    # apply update function to sliders

    S_h_diss.on_changed(update)
    S_h2.on_changed(update)
    S_h3o.on_changed(update)
    S_ho2.on_changed(update)
    S_h2o2.on_changed(update)

    # main routine #####################################################################################################
    _name = os.path.basename(argv[1])
    print("\treading file...")
    s, dt = get_system(_name)
    print("\tcomputing the difference from the non-irradiated frame...")
    shape_h2o_mols_and_find_dist(s)
    
    print("\tfinding all molecular/atomic structures...")
    species = count_h2o_mols_and_radicals(s, _h_diss_dist, _h2_creation_dist, _h3o_creation_dist,
                                          _ho2_creation_dist, _h2o2_creation_dist)

    # prepare plot

    ax.set_xlabel('time, fs')
    ax.set_ylabel('species count')

    x_axis = np.asarray(range(len(species))) * dt

    ax_h2, = ax.plot(x_axis, [len(x["h2"]) for x in species], marker='', color='y',
                     linewidth=1.5, label="$H_2$")
    ax_h202, = ax.plot(x_axis, [len(x["h2o2"]) for x in species], "-", color='m',
                       linewidth=2, label="$H_2O_2$")
    ax_h_diss, = ax.plot(x_axis, [len(x["h_diss"]) for x in species], ":", color='g',
                         linewidth=2.5, label="$H^*$")
    ax_o_diss, = ax.plot(x_axis, [len(x["o_diss"]) for x in species], "--", color='r',
                         linewidth=1, label="$O^*$")
    ax_oh_diss, = ax.plot(x_axis, [len(x["oh_diss"]) for x in species], ":", color='b',
                          linewidth=1.5, label="$OH^*$")
    ax_h3o, = ax.plot(x_axis, [len(x["h3o"]) for x in species], marker='', color='c',
                      linewidth=1.5, label="$H_3O^*$")
    ax_ho2, = ax.plot(x_axis, [len(x["ho2"]) for x in species], marker='', color='k',
                      linewidth=1, label="$HO_2^*$")

    # use 2nd axis for water count
    ax_water = ax.twinx()
    ax_h2o, = ax_water.plot(x_axis, [len(x["h2o"]) for x in species], ":", marker='', color='r',
                            linewidth=1, label="$H_2O$")
    ax_water.set_ylabel('$H_2O$ count', color="r")

    ax.legend(bbox_to_anchor=(1.13, 1))
    ax_water.legend(loc=2)

    # draw plot
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_water.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.show()
