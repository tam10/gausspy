# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:38:33 2016

@author: tam10
"""

from ase_extensions.protein_utils import Atom, Atoms, read_pdb, read_pqr, read_mol2
from pybel import readfile
from cc_utils.red_tools_utils import run_red
from cc_utils.chem_utils import get_amber_params_from_dat, get_amber_params_from_frcmod
import os
from gausspy import Gaussian
import copy
import warnings
import random
import string
import numpy as np


class Protein_Parameterisation(object):
    def __init__(self, atoms):
        
        if not hasattr(atoms, "get_pdbs"):
            raise RuntimeError("atoms must be an ASE Proteins Extension object")
        
        if os.environ.get('AMBERHOME') is None:
            raise RuntimeError("AMBERHOME not set in os.environ")
        
        alpha_nums = string.ascii_letters+string.digits
        invalid_pdbs = '; '.join([': '.join([str(a.index), a.pdb]) for a in atoms if any([s not in alpha_nums for s in a.pdb])])
        if len(invalid_pdbs) > 0:
            warnings.warn(RuntimeWarning('Invalid characters found in PDB Atom Types in atoms:\n' + invalid_pdbs + '\nThese atoms will cause errors on running Gaussian'))
        
        blank_pdbs = ', '.join([str(a.index) for a in atoms if a.pdb == ''])  
        if len(blank_pdbs) > 0:
            warnings.warn(RuntimeWarning('Blank PDB Atom Types Found in atoms:\n' + blank_pdbs))   
            
        self.Utils = self._utils(self)
        self.initial_atoms = atoms
        
        self.params = {"NSR Resnums": None,
                       "NSR Names": "CR",
                       #Protonation Options
                       "Protonation Directory": "protonation/",
                       "pH": None,
                       "Force Capping": True,
                       "Optimised Capping": True,
                       "NSR Protonation Function": self.Utils.pybel_protonate,
                       "pdb2pqr Path": "/Users/tam10/Utilities/pdb2pqr-osx-bin64-2.1.0/",
                       "pdb2pqr Options": ["--ff=amber", "--ffout=amber"],
                       #RESP Options
                       "RESP Calculation Directory": "charges/",
                       "Overwrite RESP Calculation": False,
                       #ONIOM Options
                       "Model Extraction Directory": "model/",
                       "Initial Model Atom Numbers": None,
                       "Link Atom Map": {"C": "HC",
                                         "N": "H",
                                         "S": "H"},
                       #Parameterisation Options
                       "Parameterisation Directory": "parameters/",
                       "Parameterisation Databases": [{"name": "parm99", "type": "dat"}, 
                                     {"name": "ff99SB", "type": "frcmod"}]
                      }
        self._hidden_params = {"Merged Connected NSR Resnum Groups": [],
                               "RED mol2 Filename": "Mol_m1-o1-sm.mol2",
                               "Final Model Atom Numbers": None,
                               "Original Directory": os.getcwd(),
                               "Interface Nitrogens": []}
        
        self.mm_parameters = {}
        self.nsrs = {}
        self.protonated_nsrs = {}
        self.pc_nsrs = {}
        self.protonated_atoms = None
        self.model_region = None
        
    class _utils(object):
        def __init__(self, protein_parameterisation):
            self._p = protein_parameterisation
            
        def detect_nsrs(self, pqr_file):
            pqr_file = pqr_file.split(".")[0]
            with open(pqr_file + '.pqr', 'r') as fileobj:
                filestr = fileobj.read()
            if '(omitted below):\n' in filestr and 'This is usually due to' in filestr:
                nsr_detect_lines = filestr.split('(omitted below):\n')[1].split('This is usually due to')[0].split('\n')[:-1]
                nsr_indices = list(set([int(ndl.split()[-1]) for ndl in nsr_detect_lines]))
                return sorted(nsr_indices)
            if 'WARNING: Unable to debump' in filestr:
                warn_res_lines = [l.strip() for l in filestr.split('\n') if 'WARNING: Unable to debump' in l]
                warn_resnums = [[r.split()[-1],r.split()[-3]] for r in warn_res_lines]
                warn_res_str = ", ".join(["{r} ({rn})".format(r = r, rn = rn) for r, rn in warn_resnums])
                warn_str = "NSR detection failed for residues: {r}.\nRename these residues with non-standard names to fix.".format(r = warn_res_str)
                warnings.warn(warn_str)
            else:
                return []

        def pybel_protonate(self, atoms):
            name = atoms.info.get("name")
            if name is None:
                name = atoms.get_chemical_formula()
            atoms.write_pdb(name + ".pdb")
            try:
                prot_p = readfile('pdb', name + ".pdb").next()
            except StopIteration:
                raise IOError('File {f} has no atoms')
            prot_p.addh()
            prot_p.write('pdb', name + "_pybel_protonated.pdb", overwrite = True)
            prot = read_pdb(name + "_pybel_protonated.pdb")
            prot.calculate_ambers_pdbs()
            return prot
            
        def pdb2pqr_protonate(self, atoms, filename):
            filename = filename.split()[0]
            atoms.write_pdb(filename + '.pdb')

            ph_str = '--ph-calc-method=propka --with-ph={p} '.format(p = self._p.params["pH"]) if self._p.params["pH"] is not None else ''
            os.system('{d}pdb2pqr -v {s} {p}{i} {o}'.format(d = self._p.params["pdb2pqr Path"],
                                                            s = " ".join(self._p.params["pdb2pqr Options"]), 
                                                            p = ph_str,
                                                            i = filename + '.pdb', 
                                                            o = filename +'.pqr'))
            return read_pqr(filename + '.pqr')
        
        def calculate_partial_charges(self, atoms, nsr_name):
            try:
                self._cd_in(self._p.params["RESP Calculation Directory"])
                name = atoms.info["name"]
                if name is None:
                    name = atoms.get_chemical_formula()
                rn = self._p._hidden_params["RED mol2 Filename"]
                mod_rn = nsr_name + "_partial_charges.mol2"
                rds = ['./','Data-RED/','../','../Data-RED/']
                rn_path = [rd + mod_rn for rd in rds if os.path.exists(rd + mod_rn)]
                prot_out_path = [rd + name + "-out.p2n" for rd in rds if os.path.exists(rd + name + "-out.p2n")]

                if self._p.params["Overwrite RESP Calculation"] or len(rn_path) == 0 or len(prot_out_path) == 0:
                    atoms.write_pdb(name + ".pdb")
                    os.system("Ante_RED-1.5.pl {n}".format(n = name + ".pdb"))
                    p2n_atoms = read_pdb(name + "-out.p2n")

                    ace_capping_atoms = [str(n+1) for n, a in enumerate(p2n_atoms) if a.residue == 'ACE']
                    nme_capping_atoms = [str(n+1) for n, a in enumerate(p2n_atoms) if a.residue == 'NME']

                    #Exclude ACE and NME from partial charge calculation
                    if ace_capping_atoms: ace_remark = 'REMARK INTRA-MCC 0.0 |  ' + '  '.join(ace_capping_atoms) + ' | Remove\n'
                    if nme_capping_atoms: nme_remark = 'REMARK INTRA-MCC 0.0 |  ' + '  '.join(nme_capping_atoms) + ' | Remove\n'

                    with open(name + "-out.p2n", "r") as p2n_f:
                        p2n_content = p2n_f.read()
                    split_contents = p2n_content.split('REMARK\n')

                    if ace_capping_atoms: split_contents[-2] += ace_remark
                    if nme_capping_atoms: split_contents[-2] += nme_remark
                    p2n_content = 'REMARK\n'.join(split_contents)

                    with open("mod_" + name + "-out.p2n", 'w') as p2n_f:
                        p2n_f.write(p2n_content)

                        #run_red takes ~1hr to run locally
                    os.system("sleep 1")
                    run_red("mod_" + name + "-out.p2n",nodes = 4)
                    if not os.path.exists("Data-RED/" + rn):
                        raise RuntimeError("Partial Charge calculation failed. Scratch directory must be empty")
                    os.system('cp Data-RED/{r} {m}'.format(r = rn, m = mod_rn))
                else:
                    os.system('cp {r} .'.format(r = rn_path[0]))
                    os.system('cp {p} .'.format(p = prot_out_path[0]))

                nsr_charges = read_mol2(mod_rn, atom_col_1 = 'pdb', atom_col_2 = 'amber')
                nsr_amber = read_pdb(name + "-out.p2n")
                with open(prot_out_path[0], 'r') as p2n_f:
                    pdbs = [l.strip().split()[-1] for l in p2n_f.readlines() if l.startswith('ATOM')]
                nsr_amber.calculate_ambers_pdbs()
                nsr_amber.set_pdbs(pdbs)
                
                nsr_amber.take(residues = nsr_name)
                for a in nsr_amber:
                    for b in nsr_charges:
                        if a.pdb == b.pdb:
                            a.amber_charge = b.amber_charge
                            
                return nsr_amber
            finally:
                self._cd_out()
                
        def get_oniom_model_from_indices(self, atoms, oniom_list):
            """
            Input target as an ase Atoms object
            Input indices as a list corresponding to the indices of the model region
            Returns ase Atoms object corresponding to the oniom model region
            """

            oniom_extraction = atoms.take(indices_in_tags = True)    
            indices = copy.deepcopy(oniom_list)

            oniom_extraction.set_calculator(Gaussian(label = 'model_extraction', basis = "oniom", method = "oniom(HF:UFF)=OnlyInputFiles"))
            oniom_extraction.calc.oniom_coord_params['layers'] = [indices]
            oniom_extraction.calc.oniom_coord_params['layer_mults'] = [1, 1]
            oniom_extraction.calc.set_job(time=72 ,nodes=16, memory=48000)
            oniom_extraction.calc.start(frc = True)

            atoms = oniom_extraction[oniom_list].take(indices_in_tags = True)
            with open("model_extraction.log","r") as logfile:
                lines = logfile.readlines()

            try:
                atoms_factors = [[int(w)-1 for w in l.split() if w.isdigit()] + [float(l.split()[-2])] for l in lines if "ONIOM: Cut" in l]
            except(ValueError):
                atoms_factors = [[int(w)-1 for w in l.split() if w.isdigit()] + [float(l.split()[-1].replace('\n',''))] for l in lines if "ONIOM: Cut" in l]


            for link_ind, model_ind, scale_factor in atoms_factors:
                distance = oniom_extraction.get_distance(link_ind, model_ind) * scale_factor
                oniom_extraction.set_distance(a0 = model_ind, a1 = link_ind, distance = distance, fix = 0)
                link_atom = oniom_extraction[link_ind]
                atoms += Atom(symbol = 'H',
                              position = link_atom.position,
                              tag = link_atom.tag,
                              residue = 'LNK',
                              resnum = link_atom.resnum,
                              amber_charge = link_atom.amber_charge,
                              pdb = 'H')


            self.set_link_ambers(atoms)
            return atoms

        def detect_sp2_groups(self, atoms, index=0, return_mode="atoms"):
            """Locates groups of sp2 hybridised atoms within a set of atoms and sorts them by group size, largest to smallest.
            return_mode: 'atoms', 'tags', 'indices'"""
            if isinstance(index, int):
                index = [index]
        
            name = "chrom_detection"
            atoms.write_pdb(name + '.pdb')
            pybel_nsr = readfile('pdb', name + '.pdb').next()
            hybridisation = [a.hyb for a in pybel_nsr]
        
            nns = []
            for i, n in enumerate(copy.deepcopy(atoms.get_neighbours())):
                ns = []
                if atoms[i].symbol != 'H' and hybridisation[i] == 2:
                    for a in n:
                        if atoms[a].symbol != 'H' and hybridisation[a] == 2:
                            ns.append(a)
                    ns.append(i)
                    nns.append(set(ns))
            graph = []
            while len(nns)>0:
                n0, rns = nns[0], nns[1:]
                nns = []
                for rn in rns:
                    if len(n0.intersection(rn)):
                        n0 |= rn
                    else:
                        nns.append(rn)
                n0 = list(n0)
                if return_mode == "indices":
                    graph.append([len(n0), n0])
                elif return_mode == "tags":
                    group = [t for i,t in enumerate(atoms.get_tags()) if i in n0]
                    graph.append([len(group), group])
                    
            graph = sorted(graph, reverse=True)
            graph = [g[1] for i,g in enumerate(graph) if i in index]
                
            if return_mode == "atoms":
                graph = [a for b in graph for a in b]
                return(atoms[graph])
            return graph
            
        def set_link_ambers(self, atoms):
            neighbours = atoms.get_neighbours()
            for i in range(len(atoms)):
                if atoms[i].symbol == 'H' and atoms[i].amber == '':
                    neighbour_set = neighbours[i]
                    if len(neighbour_set) != 1: 
                        raise Exception("Atom {i} ({t}) has wrong number of neighbours".format(i = i, t = atoms[i].symbol))
                    neighbour = list(neighbour_set)[0]
                    neighbour_type = atoms[neighbour].symbol
                    amber = self._p.params["Link Atom Map"].get(neighbour_type)
                    if amber is None:
                        warnings.warn("Cannot assign link atom amber for atom {i} ({t}). This may lead to errors. Add amber assignment to Link Atom Map parameter")
                        atoms[i].amber = 'H'
                    else:
                        atoms[i].amber = amber

        def merge_connected_nsr_indices(self, atoms = None):
            if atoms is None:
                atoms = self._p.initial_atoms

            ni = self._p.params["NSR Resnums"]
            if not ni:
                return []

            indices_dict = {resnum: atoms.take(resnums = resnum, indices_in_tags = True).get_tags() for resnum in ni}
            neighbours = atoms.get_neighbours()
            neighbours_dict = {neighbour_resnum: [neighbours[neighbour_index] 
                               for neighbour_index in neighbour_indices] 
                               for neighbour_resnum, neighbour_indices in indices_dict.iteritems()}
            connected_residues = [{resnum, neighbours_resnum} for neighbours_resnum, neighbours_indices in neighbours_dict.iteritems() 
                                  for resnum, indices in indices_dict.iteritems()
                                  if neighbours_resnum > resnum and len(np.intersect1d([a for b in neighbours_indices for a in b], indices)) > 0]
            
            merged_connected_residues = []
            while len(connected_residues) > 0:
                this_r, other_rs = connected_residues[0], connected_residues[1:]
                connected_residues = []
                for other_r in other_rs:
                    if len(this_r.intersection(other_r)) > 0:
                        this_r |= other_r
                    else:
                        connected_residues.append(other_r)
                merged_connected_residues.append(list(this_r))
            merged_connected_residues += [[n] for n in ni if all([n not in mcr for mcr in merged_connected_residues])]
            
            return merged_connected_residues
        
        def combine_nsrs_srs(self, original_atoms, sr, nsrs, nsr_resnum_groups, nsr_residue_names):
            if not isinstance(original_atoms, Atoms) or not hasattr(original_atoms, "get_pdbs"):
                raise RuntimeError("original_atoms must be an Atoms object")
                
            if not isinstance(sr, Atoms) or not hasattr(sr, "get_pdbs"):
                raise RuntimeError("sr must be an Atoms object")

            if not isinstance(nsrs, list):
                nsrs = [nsrs]
            for nsr in nsrs:
                if not isinstance(nsr, Atoms) or not hasattr(nsr, "get_pdbs"):
                    raise RuntimeError("nsrs must be an Atoms object or a list of Atoms objects")
            
            nrgs = nsr_resnum_groups
            if not isinstance(nrgs, list):
                nrgs = [nrgs]
            
            nrns = nsr_residue_names
            if not isinstance(nrns, list):
                nrns = [nrns]
            
            if not len(nsrs) == len(nrgs) == len(nrns):
                raise RuntimeError("Length of nsrs ({nsrs}), nsr_resnum_groups ({nrgs}) and nsr_group_names ({nrns}) must be the same".format(nsrs = len(nsrs), nrgs = len(nrgs), nrns = len(nrns)))
                
            for i in range(len(nrgs)):
                #This function corrects gaps in resnums if groups were merged in nsrs
                sr.merge_residues(nrgs[i], nrns[i])
                
                #Get tags from original atoms to align nsrs by tags with sr
                tags = list(original_atoms.take(resnums = nrgs[i]).get_tags())
                for a in nsrs[i]:
                    if a.symbol != 'H':
                        a.tag = tags.pop(0)
                        
            #Get tags from original atoms to align srs
            tags = list(original_atoms.remove(residues = [a for b in nrgs for a in b]).get_tags())
            for i in range(len(original_atoms)):
                if sr[i].symbol != 'H' and sr[i] == original_atoms[-len(tags)]:
                    try:
                        sr[i].tag = tags.pop(0)
                    except(IndexError):
                        sr[i].tag = 0
                        
            #Assemble combined protein from sr and nsr pieces
            combined = sr.take()
            for i in range(len(nrgs)):
                before_nsr = combined.take(resnums = range(0, nsrs[i][0].resnum))
                after_nsr = combined.take(resnums = range(nsrs[i][-1].resnum + 1, max(combined.get_resnums()) + 1))
                combined = before_nsr + nsrs[i].take() + after_nsr
                
            return combined
            
        def fix_interface_nitrogens(self, atoms):
            
            residues = self._p.params["NSR Names"]
            
            def get_nitrogen_type(atoms, i, neighbours_list, residue):
                if atoms[i].residue == residue and any([atoms[n].residue != residue for n in neighbours_list[i]]):
                    return 1 #Interface, NSR side
                if any([atoms[n].residue == residue for n in neighbours_list[i]]) and atoms[i].residue != residue:
                    return 2 #Interface, SR side
                else:
                    return 0 #Not on interface
                
            try:
                old_neighbours = atoms._neighbours
                atoms.calculate_neighbours()
                new_neighbours = atoms._neighbours
                interface_ns = []
                for residue in residues:
                    for i in range(len(atoms)):
                        if atoms[i].pdb == 'N':
                            n_type = get_nitrogen_type(atoms, i, new_neighbours, residue)
                            if n_type != 0 and len(new_neighbours[i]) == 3:
                                interface_ns.append([i] + new_neighbours[i])
                                
                for i_nns in interface_ns:
                    n, rest = i_nns[0], i_nns[1:]
                    for r in rest:
                        if atoms[r].symbol == 'H':
                            h = r
                            rest.remove(h)
                    if len(rest) != 2:
                        warnings.warn("Nitrogen ({n}) in residue ({r}) has more than one hydrogen neighbour while being attached to multiple residues".format(n = n, r = atoms[n].residue))
                    else:
                        self._p._hidden_param["Interface Nitrogens"].append(n)
                        c0, c1 = rest
                        mask = [1 if i == h else 0 for i in range(len(atoms))]
                        atoms.set_dihedral([c0, c1, n, h], np.pi, mask)
                        atoms.set_angle([c0, n, h], np.pi, mask)
                        atoms.set_angle([c1, n, h], 2 * np.pi / 3, mask)
            except:
                atoms._neighbours = copy.deepcopy(old_neighbours)
                
            atoms._neighbours = copy.deepcopy(old_neighbours)
            for i_nns in interface_ns:
                for i_nn in i_nns:
                    atoms._neighbours[i_nn] = new_neighbours[i_nn]

        def get_database_parameters(self, database_name, database_type):

            if database_type == 'dat':
                return get_amber_params_from_dat(database_name)
            elif database_type == 'frcmod':
                return get_amber_params_from_frcmod(database_name)
            else:
                raise RuntimeError("Database type not understood")
                
        def _cd_in(self, directory):
            self._p._hidden_params["Original Directory"] = os.getcwd()
            if not os.path.exists(directory):
                os.mkdir(directory)
            os.chdir(directory)
        
        def _cd_out(self):
            os.chdir(self._p._hidden_params["Original Directory"])
        
    def auto_protonate(self):
        
        try:    
            nns = self.params["NSR Names"]
            self.Utils._cd_in(self.params["Protonation Directory"])
            prot_nsrs = []
            
            self.prot_sr = self.Utils.pdb2pqr_protonate(self.initial_atoms, 'nsr_detection')
            self.params["NSR Resnums"] = self.Utils.detect_nsrs('nsr_detection')
            mc_nsrs = self._hidden_params["Merged Connected NSR Resnum Groups"] = self.Utils.merge_connected_nsr_indices()
            
            if isinstance(nns, str):
                if len(nns) == 3 and len(mc_nsrs) == 1:
                    if nns in self.initial_atoms.get_residues():
                        raise RuntimeError("NSR Name exists within protein")
                    self.params["NSR Names"] = nns = [nns]
                elif len(nns) == 2:
                    #Append "A", "B", "C"... to 2 letter residue name
                    nns = [nns + chr(65 + i) for i in range(len(mc_nsrs))]
                    for nn in nns:
                        if any([nn == rn for rn in self.initial_atoms.get_residues()]):
                            raise RuntimeError("NSR Name exists within protein")
                else:
                    raise RuntimeError("NSR Names must be 2 characters for multiple NSRs, or 3 letters for single NSR, or a list")
            elif len(nns) != len(mc_nsrs):
                raise RuntimeError("Length of NSR Names ({n}) must equal number of merged residues ({m})".format(n = len(nns), m = len(mc_nsrs)))
            
            self.params["NSR Names"] = nns
            
            for i, nsr_group in enumerate(mc_nsrs):
                nsr_name = nns[i]
                nsr = self.nsrs[nsr_name] = self.initial_atoms.take(resnums = nsr_group)
                
                nsr.info["name"] = nsr_name
                nsr.calculate_neighbours()
                nsr.merge_residues(nsr_group, nns[i])
                nsr.cap_sites(force = self.params["Force Capping"], optimise = self.params["Optimised Capping"])
                
                prot_nsr = self.protonated_nsrs[nsr_name] = self.params["NSR Protonation Function"](nsr)
                prot_nsr.info["name"] = nsr_name + "_prot"
                prot_nsrs.append(prot_nsr.take(residues = nns[i]))
                
            self.protonated_atoms = self.Utils.combine_nsrs_srs(self.initial_atoms, self.prot_sr, prot_nsrs, self._hidden_params["Merged Connected NSR Resnum Groups"], nns)
                    
            self.protonated_atoms.info["name"] = "master"
            self.protonated_atoms.reorder_residues()
            self.protonated_atoms.renumber_residues()
            
            chains = self.protonated_atoms.get_chains()
        
            
            self.protonated_atoms.fix_hypervalent_hydrogens()
            self.Utils.fix_interface_nitrogens(self.protonated_atoms)
            temp = self.protonated_atoms.take()
            temp.calculate_ambers_pdbs()
            for i, a in enumerate(self.protonated_atoms):
                a.chain = '' if not chains[i] or isinstance(chains[i],int) or chains[i].isdigit() else chains[i]
                if a.residue not in nns:
                    if temp[i].amber == 'DU':
                        warnings.warn("{s} ({i}) had amber type 'DU'. This can indicate an error in the structure and may negatively affect parameterisation.".format(s = a.symbol, i = i))                        
                        a.amber = a.symbol
                    else:
                        a.amber = temp[i].amber
                    
            
        finally:
            self.Utils._cd_out()
            
    def get_nsr_charges(self):
        for nsr_name in self.params["NSR Names"]:
            prot_nsr = self.protonated_nsrs[nsr_name]
            self.pc_nsrs[nsr_name] = self.Utils.calculate_partial_charges(prot_nsr, nsr_name)
            
        
        for a in self.protonated_atoms:
            for nsr_name in self.params["NSR Names"]:
                if a.residue == nsr_name:
                    for b in self.pc_nsrs[nsr_name]:
                        if a == b:
                            a.amber_charge = b.amber_charge
                            
    def get_model_region(self, atoms = None):
        
        try:
            self.Utils._cd_in(self.params["Model Extraction Directory"])
        
            if atoms == None:
                atoms = self.protonated_atoms

            initial_model_indices = self.params["Initial Model Atom Numbers"]

            if len(initial_model_indices) == 0:
                raise RuntimeError("Initial Model Atom Number parameter is empty")
            mns = list(atoms.take(tags = initial_model_indices, indices_in_tags = True).get_tags())

            neighbour_list = atoms.expand_selection(mns, mode = 'bonds', expansion = 1, inclusive = False)
            h_list = [n for n in neighbour_list if atoms[n].symbol=='H']
            mns_h = mns + h_list
            self._hidden_params["Final Model Atom Numbers"] = mns_h
            for i in range(len(atoms)):
                if i in mns_h:
                    atoms[i].atom_type = 'HETATM'
                else:
                    atoms[i].atom_type = 'ATOM'
            
            model_region = self.Utils.get_oniom_model_from_indices(atoms, sorted(mns_h))
            
            if len(model_region) == 0:
                raise RuntimeError("Model region extraction job failed")
            self.model_region = model_region
                
        finally:
            self.Utils._cd_out()
    
    def add_atom_parameters(self, atoms, name):
        try:
            self.Utils._cd_in(self.params["Parameterisation Directory"])
            if not name:
                name = atoms.info.get("name")
            if not name:
                name = "".join(random.choice(string.ascii_lowercase + string.digits) for i in range(6))
            dirname = name + "_param" + os.sep

            if not os.path.exists(dirname):
                os.mkdir(dirname)
            os.chdir(dirname)


            atoms.write_mol2('{name}.mol2'.format(name = name), atom_col_1 = 'pdb', atom_col_2 = 'amber')
            os.system('parmchk2 -a Y -f mol2 -i {name}.mol2 -o {name}.frcmod'.format(name = name))
            params = self.Utils.get_database_parameters(os.getcwd() + os.sep + '{name}.frcmod'.format(name = name), 'frcmod')

            warn_str = ''    

            for atom_type in params['types']:
                if float(atom_type['mass']) == 0:
                    warn_str += 'Atom %s has no mass!\n' % atom_type['name']
                if float(atom_type['vdwRadius']) == 0 and atom_type['name'] != 'HO':
                    warn_str += 'Atom %s has no Van der Waals radius!\n' % atom_type['name']

            for t,l in params.iteritems():
                for i in l:
                    if i.get('penalty') and float(i.get('penalty'))>0:
                        warn_str += 'Penalty score of %s in %-15s: %-2s %-2s %-2s %-2s\n' % (i['penalty'] if i.get('penalty') else '',
                                                                                             t,
                                                                                             i['t1'] if i.get('t1') else '',
                                                                                             i['t2'] if i.get('t2') else '',
                                                                                             i['t3'] if i.get('t3') else '',
                                                                                             i['t4'] if i.get('t4') else '')

            if len(warn_str) > 0:                                                            
                warnings.warn('\nErrors found in %s:\n' % (name) + warn_str)

            self.mm_parameters[name] = params

        finally:
            self.Utils._cd_out()
            
    def add_database_parameters(self):
        try:
            self.Utils._cd_in(self.params["Parameterisation Directory"])
            for db in self.params["Parameterisation Databases"]:
                self.mm_parameters[db["name"]] = self.Utils.get_database_parameters(db["name"], db["type"])
            
        finally:
            self.Utils._cd_out()
        
    def auto_parameterise(self):
        self.add_database_parameters()
        for name, pc_nsr in self.pc_nsrs.iteritems():
            self.add_atom_parameters(pc_nsr, name)
        if self.model_region is not None:
            self.add_atom_parameters(self.model_region, "model_region")
            
    def get_final_model_atomnos(self):
        return self._hidden_params["Final Model Atom Numbers"]

    def find_chromophore(self, atoms=None, expand_selection=1, include_hydrogens=False, sp2_group_index=0):
        if atoms is None:
            atoms = self.initial_atoms
        chrom = self.Utils.detect_sp2_groups(atoms, sp2_group_index, return_mode = "indices")
        chrom = [a for b in chrom for a in b]
        if expand_selection:
            ns = atoms.expand_selection(chrom, mode = 'bonds', expansion = expand_selection, inclusive = False)
            ns = [n for n in ns if atoms[n].symbol != 'H']
            chrom += ns
        if include_hydrogens:
            ns = atoms.expand_selection(chrom, mode = 'bonds', expansion = 1, inclusive = False)
            ns = [n for n in ns if atoms[n].symbol == 'H']
            chrom += ns
        self.params["Initial Model Atom Numbers"] = chrom