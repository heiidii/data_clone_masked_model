import os
import sys
import glob
import json

import pyrosetta
from pyrosetta import *
import pyrosetta.distributed.dask
from pyrosetta.rosetta.protocols.docking import *
import math
from os.path import splitext, basename
from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1
from Bio import SeqIO

pyrosetta.init('--mute all')

pdb_chains_file='masif_all_chains.txt'
outpath_all = './output_pdbs'
json_file_input='dict_pdbfiles_N2001.json'
pdbfiles = json.load(open(json_file_input, 'r'))
outpath_lowest='./lowest_pdbs_withdgres'
os.makedirs(outpath_lowest, exist_ok=True)


def read_ppi_chains_file():
    pdbs_source = []
    pdbs_source_chains = []
    with open(pdb_chains_file, 'r') as f:
            lines = f.readlines()
            pdbs_source += [t.rstrip()[:4].lower() for t in lines]
            pdbs_source_chains += [t.rstrip().lower() for t in lines]
            #print([t.rstrip().lower() for t in lines])

    pdb_chains_dict = {}
    for id, id_chain in zip(pdbs_source, pdbs_source_chains):
        #print(id, id_chain)
        pdb_chains_dict[id] = id_chain
    return pdb_chains_dict

def get_interface_analyzer_basic(partner_chain_str, scorefxn,
                           pack_separated=False) \
  -> pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover:
    interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(
    )
    interface_analyzer.fresh_instance()
    interface_analyzer.set_interface(partner_chain_str)
    interface_analyzer.set_scorefunction(scorefxn)
    interface_analyzer.set_compute_interface_energy(True)
    interface_analyzer.set_compute_interface_sc(True)
    interface_analyzer.set_calc_dSASA(True)
    interface_analyzer.set_pack_separated(pack_separated)

    return interface_analyzer

def get_prop_from_pdb(pdb, prop='dG_separated'):

    f = open(pdb, 'r')
    dat = f.readlines()
    
    for i in range(len(dat)):
        if dat[i].find(prop+' ') != -1:
            prop_val = float(dat[i].split(prop)[1].rstrip())
            return prop_val

    return None


def setup_ppi_interface_runs(overwrite=False, maxN=6000):
    pdb_chains_dict = read_ppi_chains_file()
    scorefxn = get_fa_scorefxn()
    lowest_dg_pdbs = {}
    outf=open('.done', 'w')
    for i, (key, value) in enumerate(pdbfiles.items()):
       if i > maxN:
           break
       p1 = pdb_chains_dict[key].split('_')[1]
       p2 = pdb_chains_dict[key].split('_')[2]
       
       outpath = '{}/{}_{}'.format(outpath_all, '%04d' %i, key)
       outpdbs = glob.glob('{}/*.pdb'.format(outpath))
       lowest_dg = 0.0
       lowest_pose = None
       lowest_pdb = ''
       lowest_pdb_out = '{}/{}_withbfac.pdb'.format(outpath_lowest, key)
       if os.path.exists(lowest_pdb_out):
           outf.write('{},{}\n'.format(i,key))
           continue
       if len(outpdbs)==10:
           # now we can score all; find lowest
           for cur_pdb in outpdbs:
              pose = pose_from_file(cur_pdb)
              #score = scorefxn(pose)
              #dg = pyrosetta.rosetta.core.pose.get_all_score_line_strings(pose)
              dg = get_prop_from_pdb(cur_pdb, prop='dG_separated')
              #print(dg)
              if dg < lowest_dg:
                  lowest_dg = dg
                  lowest_pose = pose.clone()
                  lowest_pdb = cur_pdb
       lowest_dg_pdbs[key] = lowest_pdb
       print(key, lowest_dg_pdbs[key])
       if lowest_pose is None:
           continue
       partners='{}_{}'.format(p1.upper(), p2.upper())
       interface_analyzer = get_interface_analyzer_basic(partners, scorefxn)
       interface_analyzer.apply(lowest_pose)
       dg_per_res = interface_analyzer.get_all_per_residue_data()
       list_dg = list(dg_per_res.dG)
       print(key, lowest_pose.size(), len(list_dg))
       assert lowest_pose.size() == len(list_dg)
       #for ires in range(pose.size()):
       #    residue = lowest_pose.residue(ires+1)
       #    for iatom in range(residue.natoms()):
       #       lowest_pose.pdb_info().bfactor(ires+1, iatom+1, 0.0)
       interface_res = list(interface_analyzer.get_interface_set())
       interface_bool = [0 if t+1 not in interface_res else 1 for t in range(pose.size())]
       for ires, dgres in enumerate(list_dg):
          #add dgres info to CA
          lowest_pose.pdb_info().bfactor(ires+1, 2, dgres)
          #interface or not info to N
          lowest_pose.pdb_info().bfactor(ires+1, 1, interface_bool[ires])
       lowest_pose.dump_file(lowest_pdb_out)
       #interface_res = list(interface_analyzer.get_interface_set())
       #dg_intres = [list_dg[i-1] for i in interface_res]

    json.dump(lowest_dg_pdbs, open('{}/lowest_dg_pdbs_N{}.json'.format(outpath_lowest, i), 'w'))
    outf.close()


setup_ppi_interface_runs()













