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
import logging
import random
import string
import numpy as np


class Protein_Parameterisation(object):
    def __init__(self, atoms):
        
        self.Utils = self._utils(self)
        
        self.params = {"NSR Resnums": None,
                       "NSR Names": "CR",
                       #Protonation Options
                       "Protonation Directory": "protonation/",
                       "pH": None,
                       "Cap Using Neighbours": True,
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
                               "Interface Nitrogens": [],
                               "Minimum Chromophore Size": 8,
                               "NSR Charges Multiplicities": None}
        

        
        self.logger = logging.getLogger(__name__)
        logging.basicConfig(format = "%(asctime)s (%(levelname)s): %(message)s", 
                            datefmt = "%d/%m/%Y %H:%M:%S", 
                            level = "INFO",
                            filename = "protein_parameterisation.log")
        
        if not hasattr(atoms, "get_pdbs"):
            raise RuntimeError("atoms must be an ASE Proteins Extension object")
        
        if os.environ.get('AMBERHOME') is None:
            raise RuntimeError("AMBERHOME not set in os.environ")
        
        alpha_nums = string.ascii_letters + string.digits
        invalid_pdbs = '; '.join([': '.join([str(a.index), a.pdb]) for a in atoms if any([s not in alpha_nums for s in a.pdb])])
        if len(invalid_pdbs) > 0:
            logging.warning('Invalid characters found in PDB Atom Types in atoms:\n' + invalid_pdbs + '\nThese atoms will cause errors on running Gaussian')
        
        blank_pdbs = ', '.join([str(a.index) for a in atoms if a.pdb == ''])  
        if len(blank_pdbs) > 0:
            logging.warning('Blank PDB Atom Types Found in atoms:\n' + blank_pdbs)  
            
        self.initial_atoms = atoms

        self.mm_parameters = {}
        self.nsrs = {}
        self.protonated_nsrs = {}
        self.pc_nsrs = {}
        self.protonated_atoms = None
        self.model_region = None
        
        logging.info("Initialised")
        
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
                logging.info("Found the following NSRs in {f}.pqr: {n}".format(n = nsr_indices, f = pqr_file))
                return sorted(nsr_indices)
            if 'WARNING: Unable to debump' in filestr:
                warn_res_lines = [l.strip() for l in filestr.split('\n') if 'WARNING: Unable to debump' in l]
                warn_resnums = [[r.split()[-1],r.split()[-3]] for r in warn_res_lines]
                warn_res_str = ", ".join(["{r} ({rn})".format(r = r, rn = rn) for r, rn in warn_resnums])
                warn_str = "NSR detection failed for residues in {f}.pqr: {r}.\nRename these residues with non-standard names to fix.".format(r = warn_res_str, f = pqr_file)
                logging.warning(warn_str)
            logging.info("Did not find any NSRs in {f}.pqr".format(f = pqr_file))
            return []

        def merge_and_cap_using_neighbours(self, atoms, resnums, name):
            res = atoms.take(resnums = resnums, indices_in_tags = True)
            res_dict = res.get_residue_dict()
            ir = [res_dict[r] for r in resnums]
            res.merge_residues(resnums, name)
            
            logging.info("Merged residues {i} ({ir}) into {f} ({fr})".format(i  = resnums,
                                                                             ir = ir,
                                                                             f  = res[0].resnum,
                                                                             fr = res[0].residue))
            tags = res.get_tags()
            
            nas0 = atoms.expand_selection(tags, inclusive = False)
            nas1 = atoms.expand_selection(nas0, inclusive = False)
            
            ris = [a for a in nas1 if a in tags]
            nas1 = [a for a in nas1 if a not in tags]
            
            for ri in ris:
                n0s = list({n0 for n0 in atoms.get_neighbours()[ri] if n0 in nas0})
                assert(len(n0s) > 0)
                if len(n0s) > 1:
                    logging.warning("Atom {a} ({n}) in residue {r} has more than one neighbour outside its own residue".format(a = atoms[ri].pdb, n = ri, r=atoms[ri].residue))
                n1s = [[n1 for n1 in atoms.get_neighbours()[n0] if n1 in nas1] for n0 in n0s]
                n1s = list({a for b in n1s for a in b})
                nas = n0s + n1s
                
                if atoms[ri].pdb == 'N':
                    if all([p in ['C', 'CA', 'O'] for p in atoms[nas].get_pdbs()]) and len(nas) == 3:
                        cap = atoms.take()[nas]
                        res += cap
                        logging.info("Adding ACE cap from {nr} to N ({an}) in {n}".format(nr = atoms[nas[0]].resnum, an = ri, n = name))
                    else:
                        res.cap_sites([i for i,t in enumerate(res.get_tags()) if i == ri])
                elif atoms[ri].pdb == 'C':
                    if all([p in ['CA', 'N'] for p in atoms[nas].get_pdbs()]) and len(nas) == 2:
                        cap = atoms.take()[nas]
                        res += cap
                        logging.info("Adding NME cap from {nr} to C ({an}) in {n}".format(nr = atoms[nas[0]].resnum, an = ri, n = name))
                    else:
                        res.cap_sites([i for i,t in enumerate(res.get_tags()) if i == ri])
                else:
                    logging.warning("Unusual capping atoms are being used for atom {a} ({an}) in {n}.".format(a = atoms[ri].pdb, an = ri, n = name))
                    
                
            res.calculate_neighbours()
            return res

        def pybel_protonate(self, atoms):
            name = atoms.info.get("name")
            if name is None:
                name = atoms.get_chemical_formula()
            atoms.write_pdb(name + ".pdb")
            i_hn = len(atoms.take('H'))
            try:
                prot_p = readfile('pdb', name + ".pdb").next()
            except StopIteration:
                raise IOError('File {f} has no atoms')
            prot_p.addh()
            prot_p.write('pdb', name + "_pybel_protonated.pdb", overwrite = True)
            prot = read_pdb(name + "_pybel_protonated.pdb")
            f_hn = len(prot.take('H'))
            logging.info("Added {h} hydrogens to {n} with Pybel".format(h = f_hn - i_hn, n = name))
            prot.calculate_ambers_pdbs()
            return prot
            
        def pdb2pqr_protonate(self, atoms, filename):
            filename = filename.split()[0]
            atoms.write_pdb(filename + '.pdb')

            ph_str = '--ph-calc-method=propka --with-ph={p} '.format(p = self._p.params["pH"]) if self._p.params["pH"] is not None else ''
            os.system('{d}pdb2pqr -v {s} {p}{f}.pdb {f}.pqr > {f}.log'.format(d = self._p.params["pdb2pqr Path"],
                                                            s = " ".join(self._p.params["pdb2pqr Options"]), 
                                                            p = ph_str,
                                                            f = filename))
            
            with open(filename + '.pqr', 'r') as fileobj:
                filestr = fileobj.read()
            if 'WARNING: Unable to debump' in filestr:
                warn_res_lines = [l.strip() for l in filestr.split('\n') if 'WARNING: Unable to debump' in l]
                warn_resnums = [[r.split()[-1],r.split()[-3]] for r in warn_res_lines]
                warn_res_str = ", ".join(["{r} ({rn})".format(r = r, rn = rn) for r, rn in warn_resnums])
                warn_str = "pdb2pqr encountered errors in residues: {r}.\n This can lead to errors later on.".format(r = warn_res_str)
                logging.warning(warn_str)
            logging.info("Protonated {f} with pdb2pqr. See {f}.log for info".format(f = filename))
            prot_atoms = read_pqr(filename + '.pqr')
            
            return prot_atoms
        
        def calculate_partial_charges(self, atoms, nsr_name):
            try:
                self._cd_in(self._p.params["RESP Calculation Directory"])
                name = atoms.info["name"]
                if name is None:
                    name = atoms.get_chemical_formula()
                rn = self._p._hidden_params["RED mol2 Filename"]
                mod_rn = nsr_name + "_partial_charges.mol2"
                rds = ['./',nsr_name + '_Data-RED/','../','../' + nsr_name + '_Data-RED/']
                rn_path = [rd + mod_rn for rd in rds if os.path.exists(rd + mod_rn)]
                prot_out_path = [rd + name + "-out.p2n" for rd in rds if os.path.exists(rd + name + "-out.p2n")]

                if self._p.params["Overwrite RESP Calculation"] or len(rn_path) == 0 or len(prot_out_path) == 0:
                    atoms.write_pdb(name + ".pdb")
                    os.system("Ante_RED-1.5.pl {n} > {l}".format(n = name + ".pdb", l = name + ".log"))
                    p2n_atoms = read_pdb(name + "-out.p2n")

                    capping_atoms = [str(n+1) for n, a in enumerate(p2n_atoms) if a.residue != nsr_name]

                    #Exclude capping atoms from partial charge calculation
                    if capping_atoms: 
                        cap_remark = 'REMARK INTRA-MCC 0.0 |  ' + '  '.join(capping_atoms) + ' | Remove\n'

                    with open(name + "-out.p2n", "r") as p2n_f:
                        split_contents = p2n_f.readlines()

                    
                    for i, split_content in enumerate(split_contents):
                        if 'CHARGE-VALUE' in split_content:
                            split_contents[i] = "REMARK CHARGE-VALUE {c}\n".format(c = self._p._hidden_params["NSR Charges Multiplicities"][nsr_name][0])
                        elif 'MULTIPLICITY-VALUE' in split_content:
                            split_contents[i] = "REMARK MULTIPLICITY-VALUE {m}\n".format(m = self._p._hidden_params["NSR Charges Multiplicities"][nsr_name][1])
                            if capping_atoms:
                                split_contents[i] += cap_remark
                                
                    p2n_content = ''.join(split_contents)

                    with open("mod_" + name + "-out.p2n", 'w') as p2n_f:
                        p2n_f.write(p2n_content)

                        #run_red takes ~1hr to run locally
                    os.system("sleep 1")
                    if self._p.params["Overwrite RESP Calculation"]:
                        logging.info("Overwriting RESP calculation for {n}...".format(n = nsr_name))
                    else:
                        logging.info("Running RESP calculation for {n}...".format(n = nsr_name))
                    run_red("mod_" + name + "-out.p2n",nodes = 4)
                    if not os.path.exists("Data-RED/" + rn):
                        raise RuntimeError("Partial Charge calculation failed for {n}. Scratch directory must be empty".format(n = nsr_name))
                    logging.info("RESP calculation complete")
                    gauss_logs = ['Data-RED/' + f for f in os.listdir('Data-RED/') if f.startswith('JOB') and f.endswith('.log')]
                    for g_log in gauss_logs:
                        with open(g_log, 'r') as g_log_obj:
                            g_log_str = g_log_obj.read()
                        if "Normal termination" not in g_log_str:
                            logging.warning("Failed job in Data-RED directory. This might mean charges are from a previous calculation")
                            
                    os.system('cp Data-RED/{r} {m}'.format(r = rn, m = mod_rn))
                    os.rename('Data-RED', nsr_name + '_Data-RED')
                else:
                    logging.info("RESP calculation already performed for {n}...".format(n = nsr_name))
                    os.system('cp {r} .'.format(r = rn_path[0]))
                    os.system('cp {p} .'.format(p = prot_out_path[0]))

                red_charges = read_mol2(mod_rn, atom_col_1 = 'pdb', atom_col_2 = 'amber')
                nsr_charges = atoms.take(indices_in_tags = True)
                
                for i, a in enumerate(nsr_charges):
                    for b in red_charges:
                        if a.pdb == b.pdb and a.residue == nsr_name:
                            a.amber_charge = b.amber_charge
                        if a.amber == 'DU':
                            logging.warning("{s} ({t}) had amber type 'DU'. This can indicate an error in the structure and may negatively affect parameterisation.".format(s = a.symbol, t = a.tag))                        
                            if a.pdb == 'N':
                                a.amber = atoms[i].amber = 'N'
                            elif a.pdb == 'CA':
                                a.amber = atoms[i].amber = 'CT'
                            elif a.pdb[0] == 'C':
                                a.amber = atoms[i].amber = 'CA'
                            elif a.pdb[0] == 'H':
                                a.amber = atoms[i].amber = 'H'
                            else:
                                a.amber = atoms[i].amber = a.symbol
                            
                
                return nsr_charges
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
                if return_mode in ["indices", "atoms"]:
                    if len(n0) >= self._p._hidden_params["Minimum Chromophore Size"]:
                        graph.append([len(n0), n0])
                elif return_mode == "tags":
                    group = [t for i,t in enumerate(atoms.get_tags()) if i in n0]
                    if len(group) >= self._p._hidden_params["Minimum Chromophore Size"]:
                        graph.append([len(group), group])
                    
            graph = sorted(graph, reverse=True)
            logging.info("{n} chromophore candidate{p} found of size{p}: {s}".format(n = len(graph),
                                                                                 s = [len(g[1]) for g in graph],
                                                                                     p = '' if len(graph) == 1 else 's'))
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
                        logging.warning("Cannot assign link atom amber for atom {i} ({t}). This may lead to errors. Add amber assignment to Link Atom Map parameter")
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
                        logging.warning("Nitrogen ({n}) in residue ({r}) has more than one hydrogen neighbour while being attached to multiple residues".format(n = n, r = atoms[n].residue))
                    else:
                        self._p._hidden_param["Interface Nitrogens"].append(n)
                        c0, c1 = rest
                        mask = [1 if i == h else 0 for i in range(len(atoms))]
                        atoms.set_dihedral([c0, c1, n, h], np.pi, mask)
                        atoms.set_angle([c0, n, h], np.pi, mask)
                        atoms.set_angle([c1, n, h], 2 * np.pi / 3, mask)
                        logging.info("Fixed interface nitrogen ({n}) in residue ({r})".format(n = n, r = atoms[n].residue))
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
                logging.info("Created the following directory: {d}".format(d = directory))
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
            
            logging.info("Grouped NSRs {i} into {f} by connectivity".format(i = self.params["NSR Resnums"], f = mc_nsrs))
            
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
            logging.info("Using the following NSR names: {n}".format(n = nns))
            
            if self._hidden_params["NSR Charges Multiplicities"] is None:
                self._hidden_params["NSR Charges Multiplicities"] = {nns[i]: [0,1] for i in range(len(nns))}
            
            for i, nsr_group in enumerate(mc_nsrs):
                nsr_name = nns[i]

                if self.params["Cap Using Neighbours"]:
                    logging.info("Capping using neighbouring residue atoms")
                    nsr = self.nsrs[nsr_name] = self.Utils.merge_and_cap_using_neighbours(self.initial_atoms, nsr_group, nsr_name)
                    
                    nsr.info["name"] = nsr_name

                else:
                    logging.info("Capping using database capping groups")
                    nsr = self.nsrs[nsr_name] = self.initial_atoms.take(resnums = nsr_group)
                    
                    nsr.info["name"] = nsr_name
                    nsr.calculate_neighbours()
                    nsr.merge_residues(nsr_group, nns[i])
                    nsr.cap_sites(force = self.params["Force Capping"], optimise = self.params["Optimised Capping"])
                
                prot_nsr = self.protonated_nsrs[nsr_name] = self.params["NSR Protonation Function"](nsr)
                prot_nsr.info["name"] = nsr_name + "_prot"
                prot_nsrs.append(prot_nsr.take(residues = nns[i]))
                
            self.protonated_atoms = self.Utils.combine_nsrs_srs(self.initial_atoms, self.prot_sr, prot_nsrs, self._hidden_params["Merged Connected NSR Resnum Groups"], nns)
                    
            logging.info("Merged standard and non-standard residues")
            
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
                        logging.warning("{s} ({i}) had amber type 'DU'. This can indicate an error in the structure and may negatively affect parameterisation.".format(s = a.symbol, i = i))                        
                        a.amber = a.symbol
                    else:
                        a.amber = temp[i].amber
                    
            logging.info("Protonation Complete")
            
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
                            if b.amber_charge == 0:
                                logging.warning("Partial charge of 0 added to atom {s} ({i}) in {n}".format(s = a.symbol, i = a.index, n = nsr_name))
                            
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
            logging.info("Model region with {n} atoms extracted using Final Atom Model Numbers".format(n = len(model_region)))
                
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
            os.system('parmchk2 -a Y -f mol2 -i {name}.mol2 -o {name}.frcmod > {name}.log'.format(name = name))
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
                logging.warning('\nErrors found in %s:\n' % (name) + warn_str)

            self.mm_parameters[name] = params

        finally:
            self.Utils._cd_out()
            
    def add_database_parameters(self):
        try:
            self.Utils._cd_in(self.params["Parameterisation Directory"])
            for db in self.params["Parameterisation Databases"]:
                self.mm_parameters[db["name"]] = self.Utils.get_database_parameters(db["name"], db["type"])
                logging.info("Added parameters for {n}".format(n = db["name"]))
            
        finally:
            self.Utils._cd_out()
        
    def auto_parameterise(self):
        self.add_database_parameters()
        for name, pc_nsr in self.pc_nsrs.iteritems():
            self.add_atom_parameters(pc_nsr, name)
            logging.info("Added parameters for {n}".format(n = name))
        if self.model_region is not None:
            self.add_atom_parameters(self.model_region, "model_region")
            logging.info("Added parameters for {n}".format(n = "model_region"))
            
    def get_final_model_atomnos(self):
        return self._hidden_params["Final Model Atom Numbers"]

    def find_chromophore(self, atoms=None, expand_selection=1, include_hydrogens=False, sp2_group_index=0):
        if atoms is None:
            atoms = self.initial_atoms
        
        c_str = "Chromophore search within {n} atoms. ".format(n = len(atoms))
            
        if expand_selection > 1:
            e_str = "Including {e} neighbours from each terminus. ".format(e = expand_selection)
        elif expand_selection == 1:
            e_str = "Including 1 neighbour from each terminus. "
        else:
            e_str = "Not including neighbouring heavy atoms. "
            
        if include_hydrogens:
            h_str = "Including hydrogens in output. "
        else:
            h_str = "Not including hydrogens in output. "
            
        s_str = "Using sp2 group index {s}. ".format(s = sp2_group_index)
            
            
        logging.info(c_str + e_str + h_str + s_str)
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

    def get_missing_param_sites(self, atoms, params = None):
        
        param_atoms = atoms.take(indices_in_tags = True)
        
        param_atoms.set_calculator(Gaussian(label = "missing_params", basis = "", method = "amber=softfirst", geom="connectivity"))
        if params:
            param_atoms.calc.set_amber_params(params)
        param_atoms.calc.start(frc = True)
        para
        
        with open("missing_params.log", "r") as mp_obj:
            mp_lines = mp_obj.readlines()
            
        allAtoms = set()
        types = []
        bondLengths  = []
        bondAngles  = []
        resnums = set()
        
        for i, mp_line in enumerate(mp_lines):
            if "Missing atomic parameters for" in mp_line:
                try:
                    a0 = int(mp_line.split()[-3]) - 1
                    
                    types.append(a0)
                    allAtoms |= {a0}
                    resnums  |= {atoms[a0].resnum}
                except(ValueError):
                    print("Could not resolve: " + mp_line)
            elif "Include all MM classes" in mp_line:
                mp_lines = mp_lines[i:]
                break
                
        for mp_line in mp_lines:
            try:
                if "Bondstretch undefined" in mp_line:
                    a0 = int(mp_line.split()[-3]) - 1
                    a1 = int(mp_line.split()[-2]) - 1
                    
                    bondLengths += [[a0, a1]]
                    allAtoms |= {a0, a1}
                    resnums  |= {atoms[a0].resnum, atoms[a1].resnum}
                    
                elif "Angle bend  undefined" in mp_line:
                    a0 = int(mp_line.split()[-4]) - 1
                    a1 = int(mp_line.split()[-3]) - 1
                    a2 = int(mp_line.split()[-2]) - 1
                    
                    bondAngles.append([a0, a1, a2])
                    allAtoms |= {a0, a1, a2}
                    resnums  |= {atoms[a0].resnum, atoms[a1].resnum, atoms[a2].resnum}
                    
            except(ValueError):
                print("Could not resolve: " + mp_line)
                
        return {"types": types, "bondLengths": bondLengths, "bondAngles": bondAngles, "resnums": list(resnums), "allAtoms": list(allAtoms)}
