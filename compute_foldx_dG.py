#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
from tqdm import tqdm
import shutil

PROJ_DIR = os.path.split(__file__)[0]
CACHE_DIR = os.path.join(PROJ_DIR, '__cache__')
if not os.path.exists(CACHE_DIR):
    os.makedirs(CACHE_DIR)
FOLDX_BIN = os.path.join(PROJ_DIR, 'foldx5', 'foldx_20251231')
if not os.path.exists(FOLDX_BIN):
    print(f'FoldX not found at {FOLDX_BIN}!')


def foldx_minimize_energy(pdb_path, out_path):
    filename = os.path.basename(os.path.splitext(pdb_path)[0]) + '_foldx.pdb'
    tmpfile = os.path.join(CACHE_DIR, filename)
    shutil.copyfile(pdb_path, tmpfile)
    p = os.popen(f'cd {CACHE_DIR}; {FOLDX_BIN} --command=Optimize --pdb={filename}')
    p.read()
    p.close()
    os.remove(tmpfile)
    filename = 'Optimized_' + filename
    tmpfile = os.path.join(CACHE_DIR, filename)
    shutil.copyfile(tmpfile, out_path)
    os.remove(tmpfile)


def foldx_dg(pdb_path, rec_chains, lig_chains):
    filename = os.path.basename(os.path.splitext(pdb_path)[0]) + '_foldx.pdb'
    tmpfile = os.path.join(CACHE_DIR, filename)
    shutil.copyfile(pdb_path, tmpfile)
    rec, lig = ''.join(rec_chains), ''.join(lig_chains)
    p = os.popen(f'cd {CACHE_DIR}; {FOLDX_BIN} --command=AnalyseComplex --pdb={filename} --analyseComplexChains={rec},{lig}')
    aff = float(p.read().split('\n')[-8].split(' ')[-1])
    p.close()
    os.remove(tmpfile)
    return aff

if __name__ == '__main__':
    # 单个样本处理
    # args = {
    #     'pdb_path': '/public/compute_dG/rank001_5EOK_619.pdb',
    #     'rec_chains': ['B'],
    #     'lig_chains': ['A']
    # }
    # dG = foldx_dg(**args)
    # print(dG)

    # 批量处理
    input_pdb_folder = '/public/PepGLAD/'
    out_file = 'result_dG/PepGLAD_foldx_dG.jsonl'
    results = []
    for pdb in tqdm(os.listdir(input_pdb_folder), desc=f'Processing foldx dG'):
        if pdb.endswith('.pdb'):
            # print(f'------------Processing {pdb}----------------')
            args = {
                    'pdb_path': os.path.join(input_pdb_folder, pdb),
                    'rec_chains': ['A'],
                    'lig_chains': ['B']
                    }
            dG = foldx_dg(**args)
            results.append((args['pdb_path'], args['rec_chains'], args['lig_chains'], dG))

    with open(out_file, 'w') as fout:
        fout.write('pdb_path\trec_chain\tlig_chain\tfoldx_dG\n')
        for i, (pdb_path, rec_chain, lig_chain, dG) in enumerate(results):
            fout.write(f'{pdb_path}\t{rec_chain}\t{lig_chain}\t{dG}\n')

    
