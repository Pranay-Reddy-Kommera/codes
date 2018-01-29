#!/usr/bin/env python

import glob
import os
import shutil

def fortran_logical_string(value):
    return '.true.' if value else '.false.'

def make_substitutions(line, opts):
    line = line.replace('__TESTCASE__', opts['testcase'])
    line = line.replace('__DIM_SPLIT__', fortran_logical_string(opts['dimsplit']))
    line = line.replace('__NOFLUX_BC__', fortran_logical_string(opts['nofluxbc']))
    line = line.replace('__NEX__', str(opts['nex']))
    line = line.replace('__NEY__', str(opts['ney']))
    line = line.replace('__NEZ__', str(opts['nez']))
    line = line.replace('__NOD__', str(opts['nod']))
    line = line.replace('__DT__', str(opts['dt']))
    line = line.replace('__SAVE_MOVIE__', fortran_logical_string(str(opts['savemovie'])))
    line = line.replace('__FRAME_FREQ__', str(opts['framefreq']))
    return line

def create_basic_mod_with_substitutions(opts):
    infile = "basic_mod.f90.template"
    outfile = "basic_mod.f90"
    with open(outfile, "wt") as fout:
        with open(infile, "rt") as fin:
            for line in fin:
                line = make_substitutions(line, opts)
                fout.write(line)

def create_rundir_and_copy_files(rundir):
    os.mkdir(rundir)
    files_to_copy = glob.glob("*.f90")
    files_to_copy.append("basic_mod.f90.template")
    files_to_copy.append("makefile")
    files_to_copy.append("submit.bsub")
    for f in files_to_copy:
        shutil.copy(f, os.path.join(rundir, f))

def prepare_rundir_and_start_run(rundir, runcmd, opts):
    print "preparing rundir = " + rundir
    create_rundir_and_copy_files(rundir)
    cwd = os.getcwd()
    os.chdir(rundir)
    create_basic_mod_with_substitutions(opts)
    os.mkdir("data") # For output data
    os.mkdir("err") # Yellowstone needs this
    os.mkdir("out") # Yellowstone needs this
    os.system(runcmd)
    os.chdir(cwd)


def get_common_options():
    opts = {} # empty dictionary
    opts['testcase'] = 'testcase_id_solid_body_rotation'
    opts['dimsplit'] = False
    opts['nofluxbc'] = True
    opts['nex'] = 16
    opts['ney'] = 16
    opts['nez'] = 16
    opts['nod'] = 3
    opts['dt'] = 'twopi/(225*should_substitute_nex_here)'
    opts['savemovie'] = False
    opts['framefreq'] = 100
    return opts

if __name__ == "__main__":
    opts = get_common_options()
    run_splits = [True, False]
    run_reslns = [8, 12]
    #runcmd = "make clean; make; ./main" # laptop
    runcmd = "make clean; make; bsub < submit.bsub" # Yellowstone
    for rs in run_splits:
      for rr in run_reslns:
        opts['dimsplit'] = rs
        opts['nex'] = rr
        opts['ney'] = rr
        opts['nez'] = rr
        opts['dt'] = 'twopi/(225*'+str(rr)+')'
        rundir = 'sbr_'+('split' if rs else '3d')+'_nod3_ne'+str(rr)
        prepare_rundir_and_start_run(rundir, runcmd, opts)

