# This file is part of the csbioscripts package
#
# Copyright (c) 2024 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import os
import json
import subprocess
from csbioscripts.fetchfromdb import downloadpage
import Bio.PDB as bpdb


class PDBentry:
    def __init__(self, pdbid):
        self.id = pdbid
        self.resolution = None
        self.expmethod = None
        self.filepath = None
        self.biopystruct = None
        self.uniprotmappings = None
        
    def fetchpdbidresolution(self, directory=None):
        jsonfile = downloadpage('https://data.rcsb.org/rest/v1/core/entry', self.id, directory=directory)
            
        if os.path.exists(jsonfile):        
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if 'rcsb_entry_info' in buf:
                    if 'resolution_combined' in buf['rcsb_entry_info']: 
                        cres = (buf['rcsb_entry_info']['resolution_combined'][0])
                        expmethod = buf['rcsb_entry_info']['experimental_method']
                        self.resolution, self.expmethod = cres, expmethod
    
    def fetchbiopythonstructure(self, pdbfile=None, directory=None):
        if pdbfile == None: # if local path not specified, fetch from RCSB PDB
            outname = f'{self.id}.cif'
            pdbfile = downloadpage('https://files.rcsb.org/download', outname, directory=directory, filename=outname)
        
        # open PDB or MMCIF parser, depending on file extension
        base, ext = os.path.splitext(pdbfile)
        parser = bpdb.MMCIFParser(QUIET=True)
        if ext.lower() == '.pdb':
            parser = bpdb.PDBParser(QUIET=True)
        
        # parse file
        if os.path.exists(pdbfile):
            self.filepath = pdbfile
            if os.path.getsize(self.filepath) > 0:
                try:
                    struct = parser.get_structure(self.id, pdbfile)
                    self.biopystruct = struct
                except:
                    self.biopystruct = None        
    
    def fetchuniprotmappings(self, directory=None):
        jsonfile = downloadpage('https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/', self.id, directory=directory)
        if os.path.exists(jsonfile):
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if self.id in buf:
                    mappings = buf[self.id]
                    self.uniprotmappings = mappings
                else:
                    self.uniprotmappings = []
    
#    def chainmapping(self, chainid):
#        if self.uniprotmappings == None:
#            self.fetchuniprotmappings()
#        if self.uniprotmappings != []:
             
    
    def printchainaspdb(self, chainid, directory=None, separate=False, suffix=None, ext='cif'):
        if self.biopystruct == None:
            self.fetchbiopythonstructure()
            
        outpdblist = []
        if separate:
            for inmodel in self.biopystruct:
                if not self.biopystruct[inmodel.id].__contains__(chainid):
                    print (f'Chain {chainid} not found in {self.id} model {inmodel.id}')
                else:
                    struct = bpdb.Structure.Structure(self.biopystruct.id)
                    model = bpdb.Model.Model(0)
                    model.add(self.biopystruct[inmodel.id][chainid])
                    struct.add(model)
                    outpdb = os.path.splitext(self.filepath)[0] + f'_{chainid}_{str(inmodel.id)}.{ext}'
                    if suffix != None: 
                        outpdb = os.path.splitext(outpdb)[0] + f'_{suffix}.{ext}'
                    if directory != None:
                        outpdb = os.path.join(directory, outpdb)
                    if ext == 'pdb':
                        printpdb(struct, outpdb)
                    else:
                        printcif(struct, outpdb)
                    outpdblist.append(outpdb)
        else:
            struct = bpdb.Structure.Structure(0)
            chainfound = False
            for inmodel in self.biopystruct:
                model = bpdb.Model.Model(inmodel.id)
                if self.biopystruct[inmodel.id].__contains__(chainid):
                    model.add(self.biopystruct[inmodel.id][chainid])
                    struct.add(model)
                    chainfound = True
            if chainfound:
                if directory != None:
                    outpdb = os.path.join(directory, outpdb)
                outpdb = os.path.splitext(self.filepath)[0] + f'_{chainid}.{ext}'
                if suffix != None: 
                    outpdb = os.path.splitext(outpdb)[0] + f'_{suffix}.{ext}'
                if ext == 'pdb':
                    printpdb(struct, outpdb)
                else:
                    printcif(struct, outpdb)
                outpdblist.append(outpdb)
        return outpdblist
    
    def chainscontainingligand(self, ligandname, directory=None):
        outchains = []
        if self.biopystruct == None:
            self.fetchbiopythonstructure(directory=directory)
        for model in self.biopystruct:
            for chain in model:
                found = False
                for res in chain:
                    if res.get_resname() == ligandname:
                        found = True
                        break
                    if found:
                        break
                if found:
                    outchains.append(chain.id)
        return sorted(list(set(outchains)))
    
    def chainlength(self, chainid, directory=None):
        if self.biopystruct == None:
            self.fetchbiopythonstructure(directory=directory)
        residues = [i for i in self.biopystruct[0][chainid]]
        return len(residues)
    
    def removealtloc(self, directory=None):
        if self.biopystruct == None:
            self.fetchbiopythonstructure(directory=directory)
        for model in self.biopystruct:
            for chain in model:
                for res in chain:
                    for atom in res:
                        altloc = atom.altloc.strip()
                        if altloc != '':
                            dellocs = [i.altloc for i in atom.disordered_get_list() if i.altloc != 'A']
                            for i in dellocs:
                                atom.disordered_remove(i)
    
    def removeheteroatoms(self, exceptions=[], directory=None):
        if self.biopystruct == None:
            self.fetchbiopythonstructure(directory=directory)
        newstruc = bpdb.Structure.Structure(self.biopystruct.id)
        for model in self.biopystruct:
            newmod = bpdb.Model.Model(model.id)
            for chain in model:
                newch = bpdb.Chain.Chain(chain.id)
                for res in chain:
                    if 'W' not in res.id[0] and 'H' not in res.id[0]:
                        newres = res.copy()
                        newch.add(newres)
                    elif 'H' in res.id[0]:
                        tag, ligand = res.id[0].split('_')
                        if ligand in exceptions:
                            newres = res.copy()
                            newch.add(newres)
                newmod.add(newch)
            newstruc.add(newmod)
        return newstruc


def printpdb(struct, path):
    io = bpdb.PDBIO()
    io.set_structure(struct)
    io.save(path)
    
def printcif(struct, path):
    io = bpdb.mmcifio.MMCIFIO()
    io.set_structure(struct)
    io.save(path)

def foldseekquery(pdbfile, db, exhaustive=False, alignment=2, cov=0.7, covmode=0, directory=None):
    fsout = f'{os.path.splitext(pdbfile)[0]}.m8'
    if directory != None:
        fsout = os.path.join(directory, fsout)
        
    # run foldseek on chain.pdb using 3did db
    if not os.path.exists(fsout):
        if exhaustive:
            result = subprocess.run(['foldseek', 'easy-search', pdbfile, db, fsout, 'tmp', '--exhaustive-search', '--alignment-type', str(alignment), '-c', str(cov), '--cov-mode', str(covmode), '--format-output','query,target,fident,qstart,qend,qlen,alnlen,evalue,lddt,prob,alntmscore', '--format-mode', '4'])
        else:
            result = subprocess.run(['foldseek', 'easy-search', pdbfile, db, fsout, 'tmp', '--alignment-type', str(alignment), '-c', str(cov), '--cov-mode', str(covmode), '--format-output','query,target,fident,qstart,qend,qlen,alnlen,evalue,lddt,prob,alntmscore', '--format-mode', '4'])

        if result.returncode != 0:
            print (f'Error: Foldseek failed. Check Foldseek installation, or if database {db} exists')
            return None
        else:
            print (f'Foldseek output file: {fsout}')
            mapping = parsefoldseekresults(fsout)
            return mapping
    else:
        mapping = parsefoldseekresults(fsout)
        return mapping

def parsefoldseekresults(fsout):
    mapping = []
    with open(fsout, 'r') as f:
        for line in f:
            query,target,fident,qstart,qend,qlen,alnlen,evalue,lddt,prob,alntmscore = line.split('\t')
                
            # record foldseek match
            if query != 'query':
                ddinum = int(target.split('_')[0])
                subnum = int(os.path.splitext(target.split('_')[1])[0])
                match = {'ddiid': ddinum,
                        'subunitnum': subnum,
                        'start': int(qstart),
                        'end': int(qend),
                        'evalue': float(evalue),
                        'lddt': float(lddt),
                        'fident': float(fident),
                        'tmscore': float(alntmscore),
                        'prob': float(prob)
                        }
                mapping.append(match)
    return mapping            

def interactingdomains(pdb1, pdb2, db, ddis, mintm=0):
    allmatches = []
    name1 = os.path.splitext(os.path.basename(pdb1))[0]
    name2 = os.path.splitext(os.path.basename(pdb2))[0]
    doms1 = foldseekquery(pdb1, db, exhaustive=True, alignment=2, cov=0.7, covmode=1)
    doms2 = foldseekquery(pdb2, db, exhaustive=True, alignment=2, cov=0.7, covmode=1)
    for i in doms1:
        if i['tmscore'] >= mintm:
            matches = [j for j in doms2 if j['ddiid'] == i['ddiid'] and j['subunitnum'] != i['subunitnum'] and j['tmscore']>=mintm]
            if len(matches) > 0:
                for j in matches:
                    ddixn = ddis[ddis['id'] == i['ddiid']]
                    pdbid = ddixn.iloc[0]['pdbid']
                    if i['subunitnum'] == 0:
                        pdbch1 = ddixn.iloc[0]['ch1']
                        pdbs1 = ddixn.iloc[0]['s1'][0]
                        pdbe1 = ddixn.iloc[0]['e1'][0]
                        pdbch2 = ddixn.iloc[0]['ch2']
                        pdbs2 = ddixn.iloc[0]['s2'][0]
                        pdbe2 = ddixn.iloc[0]['e2'][0]
                    else:
                        pdbch1 = ddixn.iloc[0]['ch2']
                        pdbs1 = ddixn.iloc[0]['s2'][0]
                        pdbe1 = ddixn.iloc[0]['e2'][0]
                        pdbch2 = ddixn.iloc[0]['ch1']
                        pdbs2 = ddixn.iloc[0]['s1'][0]
                        pdbe2 = ddixn.iloc[0]['e1'][0]
                    
                    ixn = {'protein1': name1,
                            'protein1_start': i['start'],
                            'protein1_end': i['end'],
                            'protein1_domainmatch': f'{pdbid}_{pdbch1}',
                            'protein1_domainmatch_start': pdbs1,
                            'protein1_domainmatch_end': pdbe1,
                            'protein1_domainmatch_tmscore': i['tmscore'],

                            'protein2': name2,
                            'protein2_start': j['start'],
                            'protein2_end': j['end'],
                            'protein2_domainmatch': f'{pdbid}_{pdbch2}',
                            'protein2_domainmatch_start': pdbs2,
                            'protein2_domainmatch_end': pdbe2,
                            'protein2_domainmatch_tmscore': j['tmscore']
                            }
                        
                    allmatches.append(ixn)
    return allmatches


    
def tmscorewrapper(refpdb, modpdb, fast=True):
    tmpout = 'tmp_tmscore.txt'
    tm, tm_scaled = None, None
    with open(tmpout, 'w') as f:
        if fast:
            result = subprocess.run(['/Users/kmcantos/Documents/installers/USalign', modpdb, refpdb, '-fast'], stdout=f)
        else:
            result = subprocess.run(['/Users/kmcantos/Documents/installers/USalign', modpdb, refpdb], stdout=f)
        if result.returncode == 0:
            tm = parsetmscoreoutput(tmpout)
        else:
            print ('TM scoring failed')
    subprocess.run(['rm', tmpout])
    return tm

def parsetmscoreoutput(tmout):
    tm = []
    try:
        with open(tmout, 'r') as f:
            for line in f:
                if line[0:8] == 'TM-score':
                    ln = line.split()
                    tm.append(float(ln[1]))
        return max(tm)
    except:
        return None

    
    
    

            