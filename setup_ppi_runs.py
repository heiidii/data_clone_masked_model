import os, glob
import sys
import json

xml_template_file='fastrelax_nojump.xml'
xml_template = open(xml_template_file, 'r').read()
flags_template_file='flags'
flags_template = open(flags_template_file, 'r').read()
condor_run_template_file='run.con'
condor_run_template = open(condor_run_template_file, 'r').read()
runpath='/home/saipooja/data_clone_masked_model/runs'
os.makedirs(runpath, exist_ok=True)
runpath_up='/home/saipooja/data_clone_masked_model/'
outpath_all='/home/saipooja/data_clone_masked_model/output_pdbs'
os.makedirs(runpath, exist_ok=True)
pdb_chains_file='/home/saipooja/data_clone_masked_model/masif_all_chains.txt'
pdb_files_path='/home/saipooja/data_clone_masked_model/truncated'
inputpath = '{}/inputs'.format(runpath)                

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

def setup_ppi_dg_runs(overwrite=False, maxN=6000):
    os.makedirs(inputpath, exist_ok=True)
    pdb_chains_dict = read_ppi_chains_file()
    pdbfiles = {}
    for pdbid in pdb_chains_dict:
       pdbfilename = '{}/{}_pp_trunc.pdb'.format(pdb_files_path, pdbid)
       if os.path.exists(pdbfilename):
          pdbfiles[pdbid] = pdbfilename
    for i, (key, value) in enumerate(pdbfiles.items()):
       if i > maxN:
           break
       p1 = pdb_chains_dict[key].split('_')[1]
       p2 = pdb_chains_dict[key].split('_')[2]
       p1_commas = '' + ','.join([t.upper() for t in p1])
       p2_commas = '' + ','.join([t.upper() for t in p2])
       #print(key, p1_commas, p2_commas)
       partners='{}_{}'.format(p1.upper(), p2.upper())
       
       xml_path = '{}/relax_{}.xml'.format(inputpath, key)
       outpath = '{}/{}_{}'.format(outpath_all, '%04d' %i, key)
       pdbfile_out = '{}/{}_0001.pdb'.format(outpath, pdbfiles[key].split('/')[-1].split('.')[0])
       if os.path.exists(xml_path) and (not overwrite) and \
        os.path.exists(pdbfile_out):
           continue
       open(xml_path, 'w').write(xml_template.format(placeholder_p1=p1_commas,
                           placeholder_p2=p2_commas,  
                           placeholder_partners=partners)
                            )
       flags_path = '{}/flags_{}'.format(inputpath, key)
       open(flags_path, 'w').write(flags_template.format(placeholder_pdbid=key,
                             placeholder_pdbfile=pdbfiles[key],
                             placeholder_out=outpath,
                             placeholder_xml=xml_path)
                            )
       con_path = '{}/run_{}.con'.format(runpath, key)
       open(con_path, 'w').write(condor_run_template.format(placeholder_pdbid=key,
                                                            placeholder_flags=flags_path))
       os.system('chmod 750 {}'.format(con_path))
    pdbfiles_dict_file = 'dict_pdbfiles_N{}.json'.format(i)
    json.dump(pdbfiles, open(pdbfiles_dict_file, 'w'))

    return pdbfiles_dict_file, pdbfiles


def submit_runs(json_file, N=1, pdbids=[], overwrite=False):
    pdbfiles = json.load(open(json_file, 'r'))
    os.makedirs('{}/submitted'.format(runpath_up), exist_ok=True)
    os.makedirs('{}/outerr'.format(runpath), exist_ok=True)
    for i, (key, values) in enumerate(pdbfiles.items()):
        if i > N:
            break
        flags_path = '{}/flags_{}'.format(inputpath, key)
        #print(key, i, '{}/{}_{}'.format(outpath_all, '%04d' %i, key), flags_path)
        outpath='{}/{}_{}'.format(outpath_all, '%04d' %i, key)
        os.makedirs('{}/{}_{}'.format(outpath_all, '%04d' %i, key), exist_ok=True)
        con_file = 'run_{}.con'.format(key)
        submit_file = '{}/submitted/submitted_{}_{}'.format(runpath_up, '%04d' %i, key)
        outpdbs = glob.glob('{}/*.pdb'.format(outpath))
        #print(outpath)
        if os.path.exists(submit_file) and (not overwrite) and \
            len(outpdbs)>=10:
            continue
        print(i, outpath, len(outpdbs))
        cmd = 'cd {}; condor_submit {} | tee {}'.format(runpath, con_file, submit_file)
        os.system(cmd)


def rsync_data_to_dst(dst, mode='relaxed_pdbs'):
   #currently just rsync relaxed pdbs
   if mode=='relaxed_pdbs':
      cmd_rsync = 'rsync --relative -rav saipooja@graylab.jhu.edu:{}./output_pdbs/*/*.pdb {}'.format(runpath_up, dst)
   print(cmd_rsync)
   os.system(cmd_rsync)

#jsonfile, pdbfiles = setup_ppi_dg_runs(overwrite=True, maxN=1000)
#submit_runs(jsonfile, overwrite=False, N=1000)
dst='./'
rsync_data_to_dst(dst)
