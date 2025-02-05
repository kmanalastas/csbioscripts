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
    
    def fetchbiopythonstructure(self, directory=None):
        outname = f'{self.id}.cif'
        pdbfile = downloadpage('https://files.rcsb.org/download', outname, directory=directory, filename=outname)        
        if os.path.exists(pdbfile):
            self.filepath = pdbfile
            parser = bpdb.MMCIFParser(QUIET=True)
            struct = parser.get_structure(self.id, pdbfile)
            self.biopystruct = struct
    
    def printchainaspdb(self, chainid, directory=None, separate=False):
        if self.biopystruct == None:
            self.fetchbiopythonstructure()
            
        outpdblist = []
        if separate:
            for inmodel in self.biopystruct:
                struct = bpdb.Structure.Structure(0)
                model = bpdb.Model.Model(0)
                model.add(self.biopystruct[inmodel.id][chainid])
                struct.add(model)
                outpdb = os.path.splitext(self.filepath)[0] + f'_{chainid}_{str(inmodel.id)}.pdb'
                if directory != None:
                    outpdb = os.path.join(directory, outpdb)
                printpdb(struct, outpdb)
                outpdblist.append(outpdb)
        else:
            struct = bpdb.Structure.Structure(0)
            for inmodel in self.biopystruct:
                model = bpdb.Model.Model(inmodel.id)
                model.add(self.biopystruct[inmodel.id][chainid])
                struct.add(model)
            if directory != None:
                outpdb = os.path.join(directory, outpdb)
            outpdb = os.path.splitext(self.filepath)[0] + f'_{chainid}.pdb'
            printpdb(struct, outpdb)
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


def printpdb(struct, path):
    io = bpdb.PDBIO()
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

def interactingdomains(pdb1, pdb2, db, mintm=0):
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
                    ixn = {'protein1': name1,
                            'ddiid1': i['ddiid'],
                            'subunitnum1': i['subunitnum'],
                            'p1start': i['start'],
                            'p1end': i['end'],
                            'p1tmscore': i['tmscore'],
                            'protein2': name2,
                            'ddiid2': i['ddiid'],
                            'subunitnum2': j['subunitnum'],
                            'p2start': j['start'],
                            'p2end': j['end'],
                            'p2tmscore': j['tmscore']
                            }
                    allmatches.append(ixn)
    return allmatches
    
def tmscorewrapper(refpdb, modpdb):
    tmpout = 'tmp_tmscore.txt'
    tm, tm_scaled = None, None
    with open(tmpout, 'w') as f:
        result = subprocess.run(['/Users/kmcantos/Documents/installers/USalign', modpdb, refpdb], stdout=f)
        if result.returncode == 0:
            tm = parsetmscoreoutput(tmpout)
        else:
            print ('TM scoring failed')
    subprocess.run(['rm', tmpout])
    return tm

def parsetmscoreoutput(tmout):
    tm = []
    with open(tmout, 'r') as f:
        for line in f:
            if line[0:8] == 'TM-score':
                ln = line.split()
                tm.append(float(ln[1]))
    return max(tm)

    
    
    

            