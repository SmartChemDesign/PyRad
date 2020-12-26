import configparser
import errno
import glob
import math
import os
import platform
import shutil
import subprocess
import sys

PATH = os.getcwd()
sh_script = dict()


def gen_packmol_inp(xyz, proj, cell_len, cell_cross, density, molmass):
    """
    Generates packmol input files
    :param xyz: xyz filename
    :param proj: projectile atom type
    :param cell_len: cell length in A
    :param cell_cross: cell crossection length in A
    :param density: density of compound in kg/m^3
    :param molmass: molecular mass of compound
    :return: None
    """

    # convert to float type
    cell_len = float(cell_len)
    cell_cross = float(cell_cross)
    density = float(density)
    molmass = float(molmass)

    # calculate cell parameters
    cell_vol = (cell_len * 1E-10) * ((cell_cross * 1E-10) ** 2)  # m^3
    solv_cnt = round(6.02214076E+23 * density * 1000 * cell_vol / molmass)

    print(f"%cell volume: {(cell_vol * 1e+27):.2f} nm^3")
    print(f"%molecule count: {solv_cnt}")

    # create temporary projectile xyz
    with open(f"{PATH}/packmol/PR_proj.xyz", "w+") as proj_xyz:
        proj_xyz.write(f"1\ntemporary projectile file for packmol input\n{proj} 0 0 0\n")

    # create input file for packmol
    with open(f"{PATH}/packmol/" + "/PR_packmol.inp", "w+") as inp:
        inp.write("seed -1\n")
        inp.write(f"tolerance 2.0\n")
        inp.write(f"filetype xyz\n")
        inp.write(f"output PR_initcell.xyz\n")

        inp.write(f"\nstructure PR_proj.xyz\n")
        inp.write(f"\tnumber 1\n")
        inp.write(f"\tcenter\n")
        inp.write(f"\tfixed {cell_cross / 2:.1f} {cell_cross / 2:.1f} 0. 0. 0. 0.\n")
        inp.write(f"end structure\n")

        inp.write(f"\nstructure {xyz}\n")
        inp.write(f"\tnumber {solv_cnt}\n")
        inp.write(f"\tinside box 0. 0. 0. {cell_cross:.1f} {cell_cross:.1f} {cell_len:.1f}\n")
        inp.write(f"end structure\n\n")


def clean_packmol_dir(envpath):
    """
    Clear directory with packmol
    :param envpath: path to copy resulting .xyz cell
    :return: None
    """
    # copy resulting .xyz to project dir
    try:
        os.replace("./PR_initcell.xyz", f"{envpath}/initcell.xyz")
    except OSError:
        print("!!!!!Can't copy resulting .xyz file! Check packmol.log!!!!!")
        exit()

    # clear the packmol directory of temporary .xyz and .inp files
    for i in glob.glob(f"{PATH}/packmol/*.xyz"):
        os.remove(i)
    for i in glob.glob(f"{PATH}/packmol/*.inp"):
        os.remove(i)


def gen_cell(envpath):
    """
    Runs packnmol
    :param envpath: path where to copy the resulting .xyz cell
    :return: None
    """
    # change working directory
    os.chdir(f"{PATH}/packmol/")

    # run packmol on Linux system
    if platform.system() == "Linux":
        p = subprocess.run([f"./packmol < ./PR_packmol.inp"], shell=True, check=True, stdout=subprocess.PIPE)
        with open("./packmol.log", "w") as log:
            log.write(p.stdout.decode())
    else:
        # run packmol on Linux system
        if platform.system() == "Windows":
            with open("./packmol.log", "a") as log:
                subprocess.call(["packmol.exe", "<", "PR_packmol.inp"], shell=True, stdout=log, stderr=log)
        else:
            # exit if system is unknown
            print("!!!!!Unknown system type! Can't run packmol!!!!!")
            exit()

    # clean packmol dir, copy resulting .xyz to project dir
    clean_packmol_dir(envpath)

    # go to script's root dir
    os.chdir(f"{PATH}/")


def xyz2qball(envpath, proj):
    """
    Converts .xyz file -> qb@ll .sys file
    :param envpath: path to project folder with xyz file
    :param proj: atom type of projectile
    :return: dictionary {uac: (...), pj:(x,y,z), cd:(a,b,c)}, where
        uac - set of unique atom types
        pj - coos of projectile
        cd - cell dimensions
    """
    # defining variables

    # list of atom parameters
    data = []
    # system's atoms parameters in one string
    datastr = ''
    # set of unique atom species
    atomtypes = set()

    # cell limits
    xmin = 0.0
    ymin = 0.0
    zmin = 0.0

    xmax = 0.0
    ymax = 0.0
    zmax = 0.0

    # create file and define its structure
    with open(f"{envpath}/initcell.xyz", "r") as xyzf:
        line = xyzf.readline()
        N = int(line.split()[0])  # line count
        xyzf.readline()  # skip comment

        for i in range(N):
            line = xyzf.readline()
            buffer = line.split()

            name = buffer[0]
            x = float(buffer[1])
            y = float(buffer[2])
            z = float(buffer[3])

            xmin = min(xmin, x)
            ymin = min(ymin, y)
            zmin = min(zmin, z)

            xmax = max(xmax, x)
            ymax = max(ymax, y)
            zmax = max(zmax, z)

            atomtypes.add(name)

            if name == proj:
                # projectile has fixed name: +prj
                data.append({'name': f"+prj",
                             'species': f"{name}_species",
                             'x': x,
                             'y': y,
                             'z': z})
                i -= 1
            else:
                data.append({'name': f"{name}{i + 1}",
                             'species': f"{name}_species",
                             'x': x,
                             'y': y,
                             'z': z})

        # shifting coordinates to zero            
        for i in range(N):
            if xmin < 0:
                data[i]['x'] -= xmin
            if ymin < 0:
                data[i]['y'] -= ymin
            if zmin < 0:
                data[i]['z'] -= zmin

        # shifting max coos
        if xmin < 0:
            xmax -= xmin
        if ymin < 0:
            ymax -= ymin
        if zmin < 0:
            zmax -= zmin

    # filling datastr with data from the list
    for i in range(N):
        datastr += f"atom {data[i]['name']} {data[i]['species']} \
{data[i]['x']:9.4f} {data[i]['y']:9.4f} {data[i]['z']:9.4f} angstrom\n"

    # write data to file
    with open(f"{envpath}/initcell.sys", "w+") as sysf:
        sysf.write(f"set cell {xmax:9.4f} 0 0\t0 {ymax:9.4f} 0\t0 0 {zmax:9.4f} angstrom\n")
        sysf.write(datastr)

    # TODO: rewrite the next part! This info can be obtained without file re-reading
    # get the projectile's starting coordinates and the cell params
    with open(f"{envpath}/initcell.sys", "r") as tempf:
        content = tempf.read()

    target = content.find("atom +prj")
    endtar = content[target:].find('\n') + target
    xyzcontent = content[target:endtar].split()

    x0 = float(xyzcontent[3])
    y0 = float(xyzcontent[4])
    z0 = float(xyzcontent[5])

    target = content.find("set cell")
    endtar = content[target:].find('\n') + target
    cellcontent = content[target:endtar].split()

    a = float(cellcontent[2])
    b = float(cellcontent[6])
    c = float(cellcontent[10])

    tempdict = {'uac': atomtypes, 'pj': (x0, y0, z0), 'cd': (a, b, c)}

    return tempdict


def gen_qball_scf_inp(envpath, nrowmax, ecut, ecutprec, gsp):
    """
    Generates qb@ll input for scf convergence stage and copies pseudos to work dir
    :param envpath: path to project directory
    :param nrowmax: optimized ompi thread count
    :param ecut: cutoff param for qball
    :param ecutprec: preconditioner cutoff param
    :param gsp: {pj:[x,y,z], cd:[a,b,c]} pj - coos of projectile, cd - cell dimensions
    :return: None
    """
    print("\tgenerating scf input")

    with open(f"{envpath}/scf.i", "w+") as f:
        f.write(f"set nrowmax {nrowmax}\n")
        f.write("set memory HUGE\n\n")

        f.write("set force_complex_wf ON\n\n")

        print("\tcopy pseudos to working directory")
        for i in gsp['uac']:
            f.write(f"species {i}_species {i}_HSCV_PBE-1.0.xml\n")
            src = f"{PATH}/pseudos/{i}_HSCV_PBE-1.0.xml"
            dest = f"{envpath}/{i}_HSCV_PBE-1.0.xml"
            try:
                shutil.copyfile(src, dest)
            except OSError:
                print(f"\t!!!!!{i}_HSCV_PBE-1.0.xml pseudo not found!!!!!")

        f.write("initcell.sys\n\n")

        f.write("set wf_dyn PSDA\nset xc PBE\nset vdw D3\n\n")

        f.write(f"set ecut {ecut} Rydberg\nset ecutprec {ecutprec} Rydberg\nset threshold_scf 1e-6 5\n\n")

        f.write("randomize_wf\n\n")

        f.write("run 0 20000\n\n")

        f.write("savesys scf.sys\n")
        f.write("save -states scf/scf\n\n")

        f.write("quit\n")

        # create folder for scf data
    try:
        os.mkdir(envpath + "/scf")
    except OSError as e:
        if e.errno != errno.EEXIST:
            print("!!!!!Can't create folder for /scf MOs!!!!!")
            raise


def gen_qball_move_inp(envpath, nrowmax, ecut, ecutprec, gsp, method, projmass, projeng, span, dt, savenum):
    """
    Generates qb@ll input for projectile moving stage
    :param envpath: path to project directory
    :param nrowmax: optimized ompi thread count
    :param ecut: cutoff param for qball
    :param ecutprec: preconditioner cutoff param
    :param gsp: {pj:[x,y,z], cd:[a,b,c]} pj - coos of projectile, cd - cell dimentions
    :param method: TDDFT propagation method
    :param projmass: mass of projectile in kg/m^3
    :param projeng: kin energy of projectile in keV
    :param span: how many times does projectile crosses PBC
    :param dt: TDDFT time step
    :param savenum: count of el density savepoints
    :return: None
    """
    iprojeng = int(projeng * 1000)

    print("\tgenerating projectile move stage input")
    print(f"\tE: {iprojeng} eV; span: {span}")
    with open(f"{envpath}/move_{iprojeng}_{span}.i", "w+") as f:
        f.write(f"set nrowmax {nrowmax}\n")
        f.write("set memory HUGE\n\n")

        f.write("set force_complex_wf ON\n\n")

        f.write("scf.sys\n\n")

        # VELOCITY CALCULATION______________________________________________________________________

        projvelo = math.sqrt(2 * projeng * 1.602176620898E-16 / projmass)  # m/s, translate keV -> J
        projvelo /= 2187691.2636433  # a.u. velocity
        print(f"\t\t%projectile velocity: {projvelo:.4f} a.u. velocity")

        x0, y0, z0 = gsp["pj"]  # extract init coos of proj
        a, b, c = gsp["cd"]

        # projectile ending coordinates
        x1 = x0 + span * a
        y1 = y0 + span * b
        z1 = c - z0

        # calculate velocity vector, convert to bohrs
        x = x1 - x0
        y = y1 - y0
        z = z1 - z0

        trajlen = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        print(f"\t\t%projectile trajectory length: {trajlen:.4f} angstrom")
        trajlen *= 1.889725989  # angstrom to bohrs

        # time step calculation
        flight_time = trajlen / projvelo  # a.u. time
        stepcnt = int(flight_time / dt)
        print(f"\t\t%step count: {stepcnt}")
        vx = x * 1.889725989 / flight_time  # A -> bohr
        vy = y * 1.889725989 / flight_time  # A -> bohr
        vz = z * 1.889725989 / flight_time  # A -> bohr

        f.write(f"set_velocity +prj {vx} {vy} {vz} atomicvelocity\n\n")
        # __________________________________________________________________________________________

        f.write("set wf_dyn PSDA\nset xc PBE\nset vdw D3\n\n")
        f.write(f"set ecut {ecut} Rydberg\nset ecutprec {ecutprec} Rydberg\n\n")

        f.write("load -states scf/scf\n\n")

        f.write(f"set atoms_dyn MD\nset wf_dyn {method}\n\n")

        f.write(f"set dt {dt} atomictime\nset TD_dt {dt} atomictime\n\n")

        # SAVING MANAGEMENT CONFIGURATION
        if savenum == 0:
            savenum = 1

        savestep = stepcnt // savenum

        j = 1
        stepsleft = stepcnt
        while stepsleft >= 2 * savestep:
            f.write(f"run {savestep} 1\n")
            f.write(f"save -vmd move{j}_{iprojeng}_{span}/move\n")
            try:
                os.mkdir(f"{envpath}/move{j}_{iprojeng}_{span}")
            except OSError as e:
                if e.errno != errno.EEXIST:
                    print("!!!!!Can't create density folder for move{j}_{iprojeng}_{span}!!!!!")
                    raise
            j += 1
            stepsleft -= savestep

        # rounding protection
        if stepsleft > 0:
            f.write(f"run {stepsleft} 1\n")
            f.write(f"save -vmd move{j}_{iprojeng}_{span}/move\n")
            try:
                os.mkdir(f"{envpath}/move{j}_{iprojeng}_{span}")
            except OSError as e:
                if e.errno != errno.EEXIST:
                    print("!!!!!Can't create density folder for move{j}_{iprojeng}_{span}!!!!!")
                    raise

        f.write(f"\nsave -states move_{iprojeng}_{span}/move\n")
        try:
            os.mkdir(f"{envpath}/move_{iprojeng}_{span}")
        except OSError as e:
            if e.errno != errno.EEXIST:
                print("!!!!!Can't create density folder for move_{iprojeng}_{span}!!!!!")
                raise
        f.write(f"savesys move_{iprojeng}_{span}.sys\n")

        f.write("\nquit\n")

    sh_script["move"].append(f"move_{iprojeng}_{span}.i")


def gen_qball_relax_inp(envpath, nrowmax, ecut, ecutprec, method, projeng, span, relax_dt, t, savenum):
    """
    Generates qb@ll input for host relaxation stage
    :param envpath: path to project directory
    :param nrowmax: optimized ompi thread count
    :param ecut: cutoff param for qball
    :param ecutprec: preconditioner cutoff param
    :param method: TDDFT propagation method
    :param projeng: kin energy of projectile in keV
    :param span: how many times does projectile crosses PBC
    :param relax_dt: TDDFT relaxation time step
    :param t: simulation time, femtoseconds
    :param savenum: count of RELAX el density savepoints
    :return: None
    """
    iprojeng = int(projeng * 1000)
    print("\tgenerating host relaxation stage input for corresponding moving stage")

    t *= 41.341374575751  # fs->a.u. time

    with open(f"{envpath}/relax_{iprojeng}_{span}.i", "w+") as f:
        f.write(f"set nrowmax {nrowmax}\n")
        f.write("set memory HUGE\n\n")

        f.write("set force_complex_wf ON\n\n")

        f.write(f"move_{iprojeng}_{span}.sys\n\n")

        f.write(f"set_velocity +prj 0 0 0 atomicvelocity\n\n")

        f.write(f"set wf_dyn {method}\nset atoms_dyn MD\nset xc PBE\nset vdw D3\n\n")

        f.write(f"set ecut {ecut} Rydberg\nset ecutprec {ecutprec} Rydberg\n\n")

        f.write(f"load -states move_{iprojeng}_{span}/move\n")

        f.write(f"\nset dt {relax_dt} atomictime\nset TD_dt {relax_dt} atomictime\n\n")

        # SAVING MANAGEMENT CONFIGURATION
        if savenum == 0:
            savenum = 1

        stepnum = int(t // relax_dt)
        print(f"\t\t%relax step count: {stepnum}")

        savestep = stepnum // savenum
        j = 1
        stepsleft = stepnum
        while stepsleft >= 2 * savestep:
            f.write(f"run {savestep} 1\n")
            f.write(f"save -vmd relax{j}_{iprojeng}_{span}/relax\n")
            try:
                os.mkdir(f"{envpath}/relax{j}_{iprojeng}_{span}")
            except OSError as e:
                if e.errno != errno.EEXIST:
                    print("!!!!!Can't create density folder for relax{j}_{iprojeng}_{span}!!!!!")
                    raise
            j += 1
            stepsleft -= savestep

        # rounding protection
        if stepsleft > 0:
            f.write(f"run {stepsleft} 1\n")
            f.write(f"save -vmd relax{j}_{iprojeng}_{span}/relax\n")
            try:
                os.mkdir(f"{envpath}/relax{j}_{iprojeng}_{span}")
            except OSError as e:
                if e.errno != errno.EEXIST:
                    print("!!!!!Can't create density folder for relax{j}_{iprojeng}_{span}!!!!!")
                    raise

        f.write(f"\nsave -states relax_{iprojeng}_{span}/relax\n")
        try:
            os.mkdir(f"{envpath}/relax_{iprojeng}_{span}")
        except OSError as e:
            if e.errno != errno.EEXIST:
                print("!!!!!Can't create density folder for relax_{iprojeng}_{span}!!!!!")
                raise
        f.write(f"savesys relax_{iprojeng}_{span}.sys\n")

        f.write("\nquit\n")

    sh_script["relax"].append(f"relax_{iprojeng}_{span}.i")


def gen_qball_inps(envpath, system, proj, relax, gened_sysprops):
    """
    :param envpath: path to project directory
    :param system: system parsegroup
    :param proj:  proj parsegroup
    :param relax: relax parsegroup
    :param gened_sysprops: cell, projectile and atomtypes dictionary
    :return: None
    """

    threadcnt = system.getint("ompi")
    # compute optimal nrowmax parameter
    if threadcnt < 128:
        NROWMAX = threadcnt
    elif threadcnt <= 4096:
        NROWMAX = threadcnt // 4
    else:
        NROWMAX = 2048

    # SCF STAGE
    gen_qball_scf_inp(envpath, NROWMAX, system.getint("ecut"), system.getint("ecutprec"), gened_sysprops)

    for E in proj["energy"].split(","):
        for S in proj["span"].split(","):
            # MOVE STAGE
            gen_qball_move_inp(envpath, NROWMAX, system.getint("ecut"), system.getint("ecutprec"),
                               gened_sysprops, system["propagator"].strip('"'), proj.getfloat("mass"),
                               float(E), int(S), proj.getfloat("dt"), proj.getint("savepoints"))

            # RELAX STAGE
            gen_qball_relax_inp(envpath, NROWMAX, system.getint("ecut"), system.getint("ecutprec"),
                                system["propagator"].strip('"'), float(E), int(S), relax.getfloat("dt"),
                                relax.getfloat("t"), relax.getint("savepoints"))


def gen_sh_scripts(envpath, _sys):
    """
    Generates several .sh scripts for different applications
    :param envpath: path to project dir
    :param _sys: system data from parser
    :return: None
    """
    with open(f"{envpath}/run_scf_only.sh", "w+") as f:
        f.write(f"mpirun -np {_sys.getint('ompi')} -x OMP_NUM_THREADS={_sys.getint('omp')} \
qball < scf.i > scf.o\n")

    with open(f"{envpath}/run_move_only.sh", "w+") as f:
        for elem in sh_script["move"]:
            f.write(f"mpirun -np {_sys.getint('ompi')} -x OMP_NUM_THREADS={_sys.getint('omp')} \
qball < {elem} > {elem.split('.')[0]}.o\n")

    with open(f"{envpath}/run_relax_only.sh", "w+") as f:
        for elem in sh_script["relax"]:
            f.write(f"mpirun -np {_sys.getint('ompi')} -x OMP_NUM_THREADS={_sys.getint('omp')} \
qball < {elem} > {elem.split('.')[0]}.o\n")

    with open(f"{envpath}/run_all.sh", "w+") as f:
        f.write(f"mpirun -np {_sys.getint('ompi')} -x OMP_NUM_THREADS={_sys.getint('omp')} \
qball < scf.i > scf.o\n")

        for elem in sh_script["move"]:
            f.write(f"mpirun -np {_sys.getint('ompi')} -x OMP_NUM_THREADS={_sys.getint('omp')} \
qball < {elem} > {elem.split('.')[0]}.o\n")

        for elem in sh_script["relax"]:
            f.write(f"mpirun -np {_sys.getint('ompi')} -x OMP_NUM_THREADS={_sys.getint('omp')} \
qball < {elem} > {elem.split('.')[0]}.o\n")


def create_env(_env):
    """
    Creates all needed files for computing...
    :param _env: name of config file
    :return: None
    """
    # ================================
    # PARSE CONFIG FILE
    # ================================
    print(f"===Processing {_env}...")
    cf = configparser.ConfigParser()
    cf.read(_env)

    system = cf["system"]
    host = cf["host"]
    proj = cf["projectile"]
    relax = cf["relaxation"]

    if system['propagator'].strip('"') not in ("TDEULER", "SOTD", "SORKTD", "FORKTD", "ETRS"):
        print("!!!!!CONFIG ERROR!!!!!")
        print("!!!!!propagator options are: TDEULER, SOTD, SORKTD, FORKTD, ETRS!!!!!")
        exit()

    # ================================
    # CREATE PROJECT DIRECTORY
    # ================================ 
    print("creating project directory")
    projpath = PATH + "/" + system['name'].strip('"')
    if os.path.isdir(projpath):
        print("!!!!!Directory already exists!!!!!")
        print("!!!!!Futher process will erace all data!!!!!")
        input(f"Press enter to continue, Ctrl+C to exit")
        shutil.rmtree(projpath, ignore_errors=True)

    try:
        os.mkdir(projpath)
    except OSError as e:
        if e.errno != errno.EEXIST:
            print("!!!!!Can't create project folder!!!!!")
            raise

    # ================================
    # PACKMOL ROUTINE
    # ================================
    print("creating packmol input file")
    src = PATH + "/" + host['xyz'].strip('"')
    dest = f"{PATH}/packmol/" + host['xyz'].strip('"')
    shutil.copyfile(src, dest)

    gen_packmol_inp(host["xyz"].strip('"'), proj["type"].strip('"'), host["cell_length"],
                    host["cell_section"], host["density"], host["molmass"])

    print("running packmol")
    gen_cell(projpath)

    # ================================
    # QBALL INPUT GENERATION
    # ================================
    gened_sysprops = xyz2qball(projpath, proj["type"].strip('"'))
    print("converting init cell xyz format to qb@all .sys format")

    print("generating qb@all input files")
    gen_qball_inps(projpath, system, proj, relax, gened_sysprops)

    print("generating .sh files for qball running")
    gen_sh_scripts(projpath, system)

    print("===DONE!")


def helper():
    """
    Prints the usage message and quits
    :return: None
    """
    print("Usage: PyRadCreator.py input1.ini input2.ini ...")
    print(".ini files, and the host .xyz files should be placed in the same directory as this script")
    exit()


# ================================
# MAIN
# ================================
if len(sys.argv) == 1:
    helper()
else:
    for env in sys.argv[1:]:
        sh_script.clear()
        sh_script["scf"] = "scf.i"
        sh_script["move"] = list()
        sh_script["relax"] = list()
        create_env(env)
    print("=======BYE=======")
