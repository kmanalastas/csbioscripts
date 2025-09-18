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
from csbioscripts.misc import distancematrix
from csbioscripts.seq import alignsequences
import Bio.PDB as bpdb
import numpy as np
from datetime import datetime

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
    
    def depositiondate(self, directory=None):
        jsonfile = downloadpage('https://data.rcsb.org/rest/v1/core/entry', self.id, directory=directory)
            
        if os.path.exists(jsonfile):        
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if 'pdbx_database_status' in buf:
                    if 'recvd_initial_deposition_date' in buf['pdbx_database_status']: 
                        depdate = datetime.fromisoformat(buf['pdbx_database_status']['recvd_initial_deposition_date'][:-5])
                        return depdate

    def entrytitle(self, directory=None):
        jsonfile = downloadpage('https://data.rcsb.org/rest/v1/core/entry', self.id, directory=directory)
            
        if os.path.exists(jsonfile):        
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if 'struct' in buf:
                    if 'title' in buf['struct']: 
                        title = buf['struct']['title']
                        return title
    
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
        jsonfile = downloadpage('https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/', self.id, directory=directory, filename=f'{self.id}_map.json')
        if os.path.exists(jsonfile):
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if self.id in buf:
                    mappings = buf[self.id]
                    self.uniprotmappings = mappings['UniProt']
                    self.parseuniprotmappings()
                else:
                    self.uniprotmappings = None
    
    def parseuniprotmappings(self):
        outmap = {}
        if self.uniprotmappings != None:
            for upid in self.uniprotmappings.keys():
                maplist = self.uniprotmappings[upid]['mappings']
                for mapping in maplist:
                    chid = mapping['chain_id']
                    outmap[chid] = upid
            self.uniprotmappings = outmap 
    
    def mapchaintouniprot(self, chainid):
        if self.uniprotmappings == None:
            self.fetchuniprotmappings()
        
        uniprotmatch = None
        if self.uniprotmappings != None:
            for upid in self.uniprotmappings.keys():
                for mapping in self.uniprotmappings[upid]['mappings']:
                    if mapping['chain_id'] == chainid:
                        uniprotmatch = upid
                        break
                if uniprotmatch != None:
                    break
        return uniprotmatch
                 
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

    def printsubstructurepdb(self, chainidlist, directory=None, separate=False, suffix=None, ext='cif', renamechains=False, outfile=None):
        if self.biopystruct == None:
            self.fetchbiopythonstructure()
            
        outpdblist = []
        concatchainids = ''
        for cid in chainidlist:
            concatchainids += str(cid)                        

        if separate:
            for inmodel in self.biopystruct:
                chainfound = [self.biopystruct[inmodel.id].__contains__(chainid) for chainid in chainidlist]
                if False in chainfound:
                    print (f'One chain of {chainidlist} not found')
                else:   
                    struct = bpdb.Structure.Structure(self.biopystruct.id)
                    model = bpdb.Model.Model(0)
                    for chainid in chainidlist:
                        model.add(self.biopystruct[inmodel.id][chainid])
                    struct.add(model)
                    
                    if renamechains:
                        changechainids(struct)
                    
                    # set filename
                    if outfile == None:
                        outpdb = os.path.splitext(self.filepath)[0] + f'_{concatchainids}_{str(inmodel.id)}.{ext}'
                    else:
                        outpdb = outfile + f'_{str(inmodel.id)}.{ext}'
                        
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
            chainfound = [self.biopystruct[inmodel.id].__contains__(chainid) for inmodel in self.biopystruct for chainid in chainidlist]
            if False in chainfound:
                print (f'One of the chains {chainidlist} not found')
            else:
                for inmodel in self.biopystruct:
                    model = bpdb.Model.Model(inmodel.id)
                    for chainid in chainidlist:
                        model.add(self.biopystruct[inmodel.id][chainid])
                    struct.add(model)
                
                if renamechains:
                    changechainids(struct)
                    
                # set filename
                if outfile == None:
                    outpdb = os.path.splitext(self.filepath)[0] + f'_{concatchainids}.{ext}'
                else:
                    outpdb = outfile + f'.{ext}'
                    
                if suffix != None: 
                    outpdb = os.path.splitext(outpdb)[0] + f'_{suffix}.{ext}'

                if directory != None:
                    outpdb = os.path.join(directory, outpdb)

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
        residues = self.chainsequence(chainid, directory=directory)
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
            
        self.biopystruct = newstruc
        return newstruc
    
    def matchresidues(self, selfmodid, selfcid, otherent, othermodid, othercid):
        ch1 = self.biopystruct[selfmodid][selfcid]
        ch2 = otherent.biopystruct[othermodid][othercid]
        
        newselfstruc = bpdb.Structure.Structure(0)
        newselfmod = bpdb.Model.Model(0)
        newselfch = bpdb.Chain.Chain(selfcid) 
        for res1 in ch1:
            if 'W' not in res1.id[0] and 'H' not in res1.id[0]: # if residue
                if ch2.__contains__(res1.id):
                    nres = res1.copy()
                    newselfch.add(nres)
            else:   # copy over hetatoms
                nres = res1.copy()
                newselfch.add(nres)
        newselfmod.add(newselfch)
        newselfstruc.add(newselfmod)
                
        newotherstruc = bpdb.Structure.Structure(0)
        newothermod = bpdb.Model.Model(0)
        newotherch = bpdb.Chain.Chain(othercid)                     
        for res2 in ch2:
            if 'W' not in res2.id[0] and 'H' not in res2.id[0]: # if residue
                if ch1.__contains__(res2.id):
                    nres = res2.copy()
                    newotherch.add(nres)
            else:   # copy over hetatoms
                nres = res2.copy()
                newotherch.add(nres)
        newothermod.add(newotherch)
        newotherstruc.add(newothermod)
        return newselfstruc, newotherstruc
    
    def chains_interacting(self, chainid1, chainid2, distthreshold=5.0, atomthreshold=3):
        interacting = False
        ncontacts = self.count_chain_interactions([chainid1], [chainid2], distthreshold=distthreshold)
        if ncontacts >= atomthreshold:
            interacting = True
        return interacting

    def count_chain_interactions(self, chainidlist1, chainidlist2, caonly=False):
        
        nresints = 0    # number of residue-residue interactions
        for chainid1 in chainidlist1:
            for chainid2 in chainidlist2:
                if (not self.biopystruct[0].__contains__(chainid1)) or (not self.biopystruct[0].__contains__(chainid2)):
                    return None
        
        if caonly:
            distthreshold=10.0
            atoms1 = np.array([atom.coord for chainid1 in chainidlist1 for res1 in self.biopystruct[0][chainid1] for atom in res1 if hasattr(atom, 'coord') and atom.name=='CA'])
            atoms2 = np.array([atom.coord for chainid2 in chainidlist2 for res2 in self.biopystruct[0][chainid2] for atom in res2 if hasattr(atom, 'coord') and atom.name=='CA'])
        
        else:
            distthreshold=5.0
            atoms1 = np.array([atom.coord for chainid1 in chainidlist1 for res1 in self.biopystruct[0][chainid1] for atom in res1 if hasattr(atom, 'coord')])
            atoms2 = np.array([atom.coord for chainid2 in chainidlist2 for res2 in self.biopystruct[0][chainid2] for atom in res2 if hasattr(atom, 'coord')])

#        print ('DEBUG (PDBentry.count_chain_interactions):', self.id, atoms1.shape, atoms2.shape)
        distmat = distancematrix(atoms1, atoms2)
        contactmap = np.where(distmat <= distthreshold, 1, 0)
        nresints = int(np.sum(contactmap))


#        for res1 in self.biopystruct[0][chainid1]:
#            for res2 in self.biopystruct[0][chainid2]:
#                atoms1 = np.array([atom.coord for atom in res1 if hasattr(atom, 'coord')])
#                atoms2 = np.array([atom.coord for atom in res2 if hasattr(atom, 'coord')])
#                print ('DEBUG (PDBentry.count_chain_interactions):', atoms1.shape, atoms2.shape)
#                if (len(atoms1) > 0) and (len(atoms2) > 0):
#                    distmat = distancematrix(atoms1, atoms2)
#                    contactmap = np.where(distmat <= distthreshold, 1, 0)
#                    if np.sum(contactmap) > 0:
#                        nresints += 1
        return nresints
    
    def chainsequence(self, chainid, directory=None):
        if self.biopystruct == None:
            self.fetchbiopythonstructure(directory=directory)
        
        ppb = bpdb.CaPPBuilder()
        residues = ''
        for pp in ppb.build_peptides(self.biopystruct[0][chainid]):
            residues += pp.get_sequence()
            
        return residues
    
    def chainstructuralhomologs(self, chainid, operator='relaxed_shape_match'):
        from rcsbapi.search import StructSimilarityQuery
        
        qry = StructSimilarityQuery(
            structure_search_type="entry_id",
            entry_id=self.id,
            structure_input_type="chain_id",
            chain_id=chainid,
            operator="relaxed_shape_match",
            target_search_space="polymer_entity_instance"
        )
        return qry
        
        
def pdbattributequery(attribute, operator, value):
    from rcsbapi.search import AttributeQuery
    
    qry = AttributeQuery(
        attribute=attribute,
        operator=operator,
        value = value
    )
    return qry
    
def printpdb(struct, path):
    io = bpdb.PDBIO()
    io.set_structure(struct)
    io.save(path)
    
def printcif(struct, path):
    print ('debug: biopystruct', struct)
    print ('debug: path', path)
    for mod in struct:
        print ('debug: mod', mod)
        for ch in mod:
            print ('debug: ch', ch)
    io = bpdb.mmcifio.MMCIFIO()
    io.set_structure(struct)
    io.save(path)

def changechainids(biopystruc):
    tmpids = [str(i) for i in range(101, 227)]
    chainids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    added = [str(i) for i in range(1,101)]
    chainids += added
    for model in biopystruc:
        for i, chain in enumerate(model):
            chain.id = tmpids[i]
    for model in biopystruc:
        for i, chain in enumerate(model):
            chain.id = chainids[i]

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

def find_foldseek_match_in_3did(fsmatch, ddis):
    ddixn = ddis[ddis['id'] == fsmatch['ddiid']]
    pdbid = ddixn.iloc[0]['pdbid']
    matchpdb = PDBentry(pdbid)
    matchpdb.fetchuniprotmappings()
    
    if fsmatch['subunitnum'] == 0:
        pdbchain = ddixn.iloc[0]['ch1']
        pdbstart = ddixn.iloc[0]['s1'][0]
        pdbend = ddixn.iloc[0]['e1'][0]
    else:
        pdbchain = ddixn.iloc[0]['ch2']
        pdbstart = ddixn.iloc[0]['s2'][0]
        pdbend = ddixn.iloc[0]['e2'][0]

#    mup, mname = matchpdb.uniprotmappings[pdbchain][0]
    mup = matchpdb.uniprotmappings[pdbchain]
    
    details = {'pdbid': pdbid,
                'pdbchain': pdbchain,
                'pdbstart': pdbstart,
                'pdbend': pdbend,
                'uniprotid': mup,
                'uniprotname': mup #mname
    }
    return details

def domainhits(pdb, db, ddis, mintm=0):
    allmatches = []
    name = os.path.splitext(os.path.basename(pdb))[0]
    doms = foldseekquery(pdb, db, exhaustive=True, alignment=2, cov=0.7, covmode=1)
    for i in doms:
        if i['tmscore'] >= mintm:
            details = {'querystart': i['start'],
                        'queryend': i['end']
            }
            match = find_foldseek_match_in_3did(i, ddis)
            details.update(match)
            details.update({'tmscore': i['tmscore'], 'fident': i['fident']})
            allmatches.append(details)
    return allmatches

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
                detailsi = find_foldseek_match_in_3did(i, ddis)
                for j in matches:
                    detailsj = find_foldseek_match_in_3did(j, ddis)
                    
                    ixn = {'protein1': name1,
                            'protein1_start': i['start'],
                            'protein1_end': i['end'],
                            'protein1_domainmatch': f"{detailsi['pdbid']}_{detailsi['pdbchain']}",
                            'protein1_domainmatch_start': detailsi['pdbstart'],
                            'protein1_domainmatch_end': detailsi['pdbend'],
                            'protein1_domainmatch_uniprotid': detailsi['uniprotid'],
                            'protein1_domainmatch_uniprotname': detailsi['uniprotname'],
                            'protein1_domainmatch_tmscore': i['tmscore'],

                            'protein2': name2,
                            'protein2_start': j['start'],
                            'protein2_end': j['end'],
                            'protein2_domainmatch': f"{detailsj['pdbid']}_{detailsj['pdbchain']}",
                            'protein2_domainmatch_start': detailsj['pdbstart'],
                            'protein2_domainmatch_end': detailsj['pdbend'],
                            'protein2_domainmatch_uniprotid': detailsj['uniprotid'],
                            'protein2_domainmatch_uniprotname': detailsj['uniprotname'],
                            'protein2_domainmatch_tmscore': j['tmscore']
                            }
                        
                    allmatches.append(ixn)
    return allmatches

def dockq(refpdb, modpdb):
    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
    
    model = load_PDB(modpdb)
    native = load_PDB(refpdb)

    # native:model chain map dictionary for two interfaces
    chain_map = {"A":"A", "B":"B"}
    # returns a dictionary containing the results and the total DockQ score
    return run_on_all_native_interfaces(model, native, chain_map=chain_map)
    
def dockqscorewrapper(refpdb, modpdb):
    # print pdb files with matching residues only
    ref = PDBentry('ref')
    mod = PDBentry('mod')
    ref.fetchbiopythonstructure(refpdb)
    mod.fetchbiopythonstructure(modpdb)
    commonchainids = [chain.id for chain in ref.biopystruct[0] if mod.biopystruct[0].__contains__(chain.id)]
    
    ## build new structures 
    refstruct = bpdb.Structure.Structure(ref.biopystruct.id)
    refmodel = bpdb.Model.Model(0)
    modstruct = bpdb.Structure.Structure(mod.biopystruct.id)
    modmodel = bpdb.Model.Model(0)
    for chainid in commonchainids:
        refchain = bpdb.Chain.Chain(chainid)
        modchain = bpdb.Chain.Chain(chainid)
        refseq = ref.chainsequence(chainid)
        modseq = mod.chainsequence(chainid)
        alignment = alignsequences(refseq, modseq)
        mapping = equivalentresidues(ref.biopystruct[0][chainid], mod.biopystruct[0][chainid], alignment[0])
        #print (chainid, mapping)
        refres = [res for res in ref.biopystruct[0][chainid] if res.id[1] in mapping]
        modres = [mod.biopystruct[0][chainid][mapping[res.id[1]]] for res in refres]
        for i in range(len(refres)):
            # change refres residue number to match modres
            newref = refres[i].copy()
            newref.id = modres[i].id
            
            # add residues to respective chains
            refchain.add(newref)
            modchain.add(modres[i])
        refmodel.add(refchain)
        modmodel.add(modchain)
    refstruct.add(refmodel)
    modstruct.add(modmodel)
    tmpref = 'tmpref.pdb'
    tmpmod = 'tmpmod.pdb'
    printpdb(refstruct, 'tmpref.pdb')
    printpdb(modstruct, 'tmpmod.pdb')
    
    # call dockq installation
    tmpout = 'tmp_dockq.txt'
    os.system(f'/Users/kmcantos/Documents/installers/DockQ-master/DockQ.py {tmpmod} {tmpref} > {tmpout}')
    dqscore = -1
    with open(tmpout, 'r') as f:
        for line in f:
            if line[0:5] == 'DockQ':
                ln = line.strip().split()
                dqscore = ln[1]
                break
    name = os.path.splitext(os.path.basename(modpdb))[0]
    print (name, dqscore)    

    # clean up intermediate files
    #os.system(f'rm {tmpref} {tmpmod} {tmpout}')

    return dqscore

def equivalentresidues(refchain, modchain, alignment):
    mapping = {}
    refresidues = [res for res in refchain]
    modresidues = [res for res in modchain]
    refind = -1
    modind = -1
    #print (alignment[0][0], len(alignment[0][0]))
    for i in range(len(alignment[0])):
#        print (i)
#        print (alignment[0][i], alignment[1][i])
        if alignment[0][i] != '-':
            refind += 1
        if alignment[1][i] != '-':
            modind += 1
        if alignment[0][i] != '-' and alignment[1][i] != '-':
            mapping[refresidues[refind].id[1]] = modresidues[modind].id[1] 
    return mapping

def tmscorewrapper(refpdb, modpdb, fast=True, multimer=False):
    tmpout = 'tmp_tmscore.txt'
    tm, tm_scaled = None, None
    with open(tmpout, 'w') as f:
        mm = '0'
        ter = '2'
        if multimer:
            mm = '1'
            ter = '0'
        if fast:
            result = subprocess.run(['/Users/kmcantos/Documents/installers/USalign', modpdb, refpdb, '-fast', '-mm', mm, '-ter', ter], stdout=f)
        else:
            result = subprocess.run(['/Users/kmcantos/Documents/installers/USalign', modpdb, refpdb, '-mm', mm, '-ter', ter], stdout=f)
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

    
    
    

            
            