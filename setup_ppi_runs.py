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




def setup_ppi_dg_runs(overwrite=False, maxN=6000, pdbfiles=None):
    os.makedirs(inputpath, exist_ok=True)
    pdb_chains_dict = read_ppi_chains_file()
    if pdbfiles is None:
       pdbfiles = {}
       for pdbid in pdb_chains_dict:
          pdbfilename = '{}/{}_pp_trunc.pdb'.format(pdb_files_path, pdbid)
          #print(pdbid, pdbfilename)
          if os.path.exists(pdbfilename):
             pdbfiles[pdbid] = pdbfilename
    #print(len(pdbfiles.keys()))
    for i, (key, value) in enumerate(pdbfiles.items()):
       if i > maxN:
           break
       p1 = pdb_chains_dict[key].split('_')[1]
       p2 = pdb_chains_dict[key].split('_')[2]
       p1_commas = '' + ','.join([t.upper() for t in p1])
       p2_commas = '' + ','.join([t.upper() for t in p2])
       #print(key, p1_commas, p2_commas)
       partners='{}_{}'.format(p1.upper(), p2.upper())
       chains='{}{}'.format(p1.upper(), p2.upper())
       cleanup_name = '{}_{}.pdb'.format(value.split('/')[-1].split('.pdb')[0], chains)
       cleanupfile_exists = os.path.exists(cleanup_name)
       xml_path = '{}/relax_{}.xml'.format(inputpath, key)
       outpath = '{}/{}_{}'.format(outpath_all, '%04d' %i, key)
       pdbfile_out = '{}/{}_0001.pdb'.format(outpath, pdbfiles[key].split('/')[-1].split('.')[0])
       if os.path.exists(xml_path) and (not overwrite) and \
        os.path.exists(pdbfile_out) and (not cleanupfile_exists):
           continue
       print('cleanup: ',key, cleanup_name)
       open(xml_path, 'w').write(xml_template.format(placeholder_p1=p1_commas,
                           placeholder_p2=p2_commas,  
                           placeholder_partners=partners)
                            )
       flags_path = '{}/flags_{}'.format(inputpath, key)
       if not cleanupfile_exists:
          open(flags_path, 'w').write(flags_template.format(placeholder_pdbid=key,
                             placeholder_pdbfile=pdbfiles[key],
                             placeholder_out=outpath,
                             placeholder_xml=xml_path)
                            )
       else:
          open(flags_path, 'w').write(flags_template.format(placeholder_pdbid=key,
                             placeholder_pdbfile=os.path.abspath(cleanup_name),
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


def submit_runs(json_file, N=1, pdbids=[], overwrite=False,
                submit=True, selected_ids=None):
    pdbfiles = json.load(open(json_file, 'r'))
    os.makedirs('{}/submitted'.format(runpath_up), exist_ok=True)
    os.makedirs('{}/outerr'.format(runpath), exist_ok=True)
    if selected_ids is None:
       selected_ids = list(pdbfiles.keys())
    cleanup_ids=[]
    for i, (key, values) in enumerate(pdbfiles.items()):
        if i > N:
            break
        if key not in selected_ids:
           continue
        flags_path = '{}/flags_{}'.format(inputpath, key)
        #print(key, i, '{}/{}_{}'.format(outpath_all, '%04d' %i, key), flags_path)
        outpath='{}/{}_{}'.format(outpath_all, '%04d' %i, key)
        os.makedirs('{}/{}_{}'.format(outpath_all, '%04d' %i, key), exist_ok=True)
        con_file = 'run_{}.con'.format(key)
        submit_file = '{}/submitted/submitted_{}_{}'.format(runpath_up, '%04d' %i, key)
        outpdbs = glob.glob('{}/*.pdb'.format(outpath))
        #print(outpath)
        if os.path.exists(submit_file) and (not overwrite) and \
            len(outpdbs)>=5:
            continue
        if submit:
         print(i, outpath, len(outpdbs))
         cmd = 'cd {}; condor_submit {} | tee {}'.format(runpath, con_file, submit_file)
         os.system(cmd)
        else:
         cleanup_ids.append('{}_{}'.format('%04d' %i, key))
    return cleanup_ids


def rsync_data_to_dst(dst, mode='relaxed_pdbs'):
   #currently just rsync relaxed pdbs
   if mode=='relaxed_pdbs':
      cmd_rsync = 'rsync --relative -rav saipooja@graylab.jhu.edu:{}./output_pdbs/*/*.pdb {}'.format(runpath_up, dst)
   print(cmd_rsync)

def remove_bad_ids():
   with open('redo_ids_paths.txt', 'r') as f:
      bad_ids = [t.rstrip() for t in f.readlines()]
   for bid in bad_ids:
      fullpath = '{}/{}'.format(outpath_all, bid)
      cmd_rm = 'rm -rf {}'.format(fullpath)
      if os.path.exists(fullpath):
        print(cmd_rm)
        os.system(cmd_rm)  

#jsonfile='dict_pdbfiles_N5397.json'
#pdbfiles_from_dict = json.load(open(jsonfile, 'r'))
#jsonfile, pdbfiles = setup_ppi_dg_runs(overwrite=False, maxN=5500, pdbfiles=pdbfiles_from_dict)
jsonfile='dict_pdbfiles_N5397.json'
cleaned_pdbs = glob.glob('./*.pdb')
#print(cleaned_pdbs)
selected_ids = []
selected_ids = [os.path.basename(t).rstrip()[:4] for t in cleaned_pdbs]
#print(selected_ids)
cleanup_ids = submit_runs(jsonfile, overwrite=False, N=5500, submit=False)
print(cleanup_ids)
print(len(cleanup_ids))
#dst='/Users/saipooja/Documents/Repositories/data_clone_masked_model'
#rsync_data_to_dst(dst)
#remove_bad_ids()
