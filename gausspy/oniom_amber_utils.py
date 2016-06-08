# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:38:33 2016

@author: tam10
"""

from ase_extensions.protein_utils import Atom, Atoms, read_pdb, read_pqr, read_mol2
from pybel import readfile
from cc_utils.red_tools_utils import run_red
from cc_utils.chem_utils import get_amber_params_from_dat, get_amber_params_from_frcmod
import xml.etree.ElementTree as tree
from xml.dom import minidom
import os
from gausspy import Gaussian
import copy
import warnings
import random
import string
import numpy as np

def protonate(atoms, 
              return_params=False, 
              non_standard_indices=None,
              non_standard_residue_name='CHR', 
              nsr_prot_function=None, 
              use_red=True, 
              ph=None, 
              force_capping=False,
              debug_level=1, 
              logfile='',
              pdb2pqr_path="/Users/tam10/Utilities/pdb2pqr-osx-bin64-2.1.0/"):
    """Protonate an atoms object with standard and non-standard residues.
    
    debug_levels:
        0:              Suppress warnings
        1:              (Default) Warnings not suppressed
        2:              Intermediate files kept
        3:              Progress printed to stdout
        4:              Returns intermediate objects in dictionary"""
    
    #Check initial atoms object
    alpha_nums = string.ascii_letters+string.digits
    invalid_pdbs = '; '.join([': '.join([str(a.index), a.pdb]) for a in atoms if any([s not in alpha_nums for s in a.pdb])])
    if len(invalid_pdbs) > 0:
        warnings.warn(RuntimeWarning('Invalid characters found in pdbs in atoms:\n'+invalid_pdbs+'\nThese atoms will cause errors on running Gaussian'))
    
    blank_pdbs = ', '.join([str(a.index) for a in atoms if a.pdb == ''])  
    if len(blank_pdbs) > 0:
        warnings.warn(RuntimeWarning('Blank pdbs found in atoms:\n'+blank_pdbs))
    
    pwd = os.getcwd()
    d = debug_level
    ni = sorted(non_standard_indices)
    nn = non_standard_residue_name
    
    if d > 4:
        r_dict = {'all' : {}, 'nsr' : {}, 'sr' : {}}
        r_dict['all']['original'] = atoms.take()
    
    rn = 'Mol_m1-o1-sm.mol2'
    dirname = 'protonation/'  
    
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        if d > 2:
            print("Created directory: {d}".format(d = dirname))
    
    os.chdir(dirname)    
    
    if not ni: ni = []
    sr = atoms.take(indices_in_tags = True)
    sr.set_atom_types(['ATOM'] * len(sr))
    if d>4:
        r_dict['sr']['original'] = sr.remove(resnums = ni)
    
    if len(ni) > 0:
        #Take nsr from original file
        nsr = sr.take(resnums = ni)
        if d > 4:
            r_dict['nsr']['original'] = nsr.take()
        if d > 2:
            print("Extracted non-standard residues: {n}".format(n = ni))
        
        #Merge nsrs and rename to non_standard_residue_name
        sr.merge_residues(ni, nn)
        
        #Merge and cap nsr
        nsr.merge_residues(ni,nn)
        if d > 4:
            r_dict['nsr']['merged'] = nsr.take()
        if d > 0:
            if len(nsr.get_ace_cap_sites()) == 0 and not force_capping:
                warnings.warn(RuntimeWarning("No ACE capping sites found. Use 'force_capping=True"))
            if len(nsr.get_nme_cap_sites()) == 0 and not force_capping:
                warnings.warn(RuntimeWarning("No NME capping sites found. Use 'force_capping=True"))
        nsr.cap_sites(force = force_capping)
        if d > 4:
            r_dict['nsr']['capped'] = nsr.take()
        if d > 2:
            print("Merged and capped non-standard residues and renamed to {n}".format(n = nn))
            
        #Protonate nsr
        nsr.write_pdb('nsr_capped.pdb')
        if nsr_prot_function is None:
            nsr_prot_function=pybel_protonate
        os.system('cp ' + nsr_prot_function('nsr_capped.pdb', (d > 2)) + ' ./nsr_capped_prot.pdb')
        os.system('Ante_RED-1.5.pl nsr_capped_prot.pdb')
        nsr_capped = read_pdb('nsr_capped_prot-out.p2n')
        
        ace_capping_atoms = [str(n+1) for n, a in enumerate(nsr_capped) if a.residue == 'ACE']
        nme_capping_atoms = [str(n+1) for n, a in enumerate(nsr_capped) if a.residue == 'NME']
        
        #Exclude ACE and NME from partial charge calculation
        ace_remark = 'REMARK INTRA-MCC 0.0 |  ' + '  '.join(ace_capping_atoms) + ' | Remove\n'
        nme_remark = 'REMARK INTRA-MCC 0.0 |  ' + '  '.join(nme_capping_atoms) + ' | Remove\n'

        with open('nsr_capped_prot-out.p2n') as p2n_f:
            p2n_content = p2n_f.read()
        split_contents = p2n_content.split('REMARK\n')
        
        split_contents[-2] += ace_remark
        split_contents[-2] += nme_remark
        p2n_content = 'REMARK\n'.join(split_contents)
        
        with open('mod_nsr_capped_prot-out.p2n', 'w') as p2n_f:
            p2n_f.write(p2n_content)
        
        if use_red:
            #run_red takes ~1hr to run locally
            if d > 2:
                print("Calculating non-standard residue RESP charges")
            os.system("sleep 1")
            run_red('mod_nsr_capped_prot-out.p2n',nodes = 4)
            os.system('cp Data-RED/{r} .'.format(r = rn))
        else:
            try:
                rds = ['./','Data-RED/','../','../Data-RED/']
                rn_path = [rd for rd in rds if os.path.exists(rd + rn)][0]
                os.system('cp {p} .'.format(p = rn_path + rn))
            except(IndexError):
                os.chdir(pwd)
                raise IOError("RED Output {r} not found. Run with 'use_red=True'".format(r = rn))
            if d > 2:
                print("Found RED output {r} in {p}".format(r = rn,p = rn_path))
        try:
            nsr_charges = read_mol2(rn, atom_col_1 = 'pdb', atom_col_2 = 'amber')
        except(IOError):
            os.chdir(pwd)
            raise IOError("RED Output {r} not found. Scratch directory must be empty".format(r = rn))
        nsr_amber = read_pdb('nsr_capped_prot-out.pdb').take(residues = nn)
        nsr_amber.calculate_ambers_pdbs(debug = (d > 1))
        nsr_amber.set_amber_charges(nsr_charges.get_amber_charges())
        nsr_amber.set_resnums(ni[0])
        
        if d > 4:
            r_dict['nsr']['final'] = nsr_amber.take()
        
        tags = list(sr.take(residues='CHR').get_tags())
        tags.reverse()
        for a in nsr_amber:
            if a.symbol != 'H':
                a.tag = tags.pop()
    
    sr.write_pdb('sr.pdb')
    tags = list(sr.get_tags())
    tags.reverse()
    ph_str = '--ph-calc-method=propka --with-ph={p} '.format(p = ph) if ph else ''
    os.system('{f}/pdb2pqr -v --ff=amber --ffout=amber --nodebump {p}sr.pdb sr.pqr > pdb2pqr.log'.format(p = ph_str, f = pdb2pqr_path))
    master = read_pqr('sr.pqr').remove(residues='CHR')
    
    if len(ni) > 0:
        master = master.take(resnums = range(0, ni[0])) + nsr_amber + master.take(resnums = range(ni[0]+1, max(master.get_resnums()) + 1))
        
    for a in master:
        if a.symbol != 'H':
            try:
                a.tag = tags.pop()
            except(IndexError):
                a.tag = 0
        
    if d > 4:
        r_dict['all']['merged'] = master
    master.info["name"] = "master"
    master.clean(debug = (d > 1))
    
    os.chdir(pwd)    
    if d < 2:
        os.system('cp {p} .'.format(p = dirname + rn))
        os.system("rm -r {d}".format(d = dirname))
    
    if d > 4:
        r_dict['all']['final'] = master
        return r_dict
    else:
        return master
        
def pybel_protonate(filename, debug = False):
    nsr_capped_prot_p = readfile('pdb', filename).next()
    nsr_capped_prot_p.addh()
    nsr_capped_prot_p.write('pdb', 'nsr_capped.pdb', overwrite = True)
    if debug:
        print("Protonated non-standard residues with Pybel")
    
    return 'nsr_capped.pdb'
    #return filename.replace('.pdb','.pdb')

def pdb2pqr_protonate(protein, pp_nsr, pdb2pqr_path, ph = None, debug = False):
    
    #Tests:
    
    #Look for repeats in pdbs which will break at AA.xml
    s = np.sort(pp_nsr.get_pdbs())
    duplicates = s[s[1:] == s[:-1]]
    if len(duplicates)>0:
        raise RuntimeError('Duplicates found in nsr pdbs: {p}'.format(p = s[s[1:] == s[:-1]]))
    
    os.system('rm -r pdb2pqr; mkdir pdb2pqr; cd pdb2pqr; cp -r {p}/* .'.format(p = pdb2pqr_path))
    
    #Check pdbs align in protein
    
    #Handle AMBER.DAT:
    
    nsr_params = get_atom_parameters(pp_nsr,'nsr')['types']
    vdws = {p['element']:p['vdwRadius'] for p in nsr_params}
    pdb2pqr_additional_params="\n".join("{r}\t{p}\t{c:.6f}\t{v:.4f}\t{m}".format(r = a.residue, 
                                                                                 p = a.pdb,
                                                                                 c = a.amber_charge,
                                                                                 v = float(vdws[a.amber]), 
                                                                                 m = a.amber) for a in pp_nsr)
    
    with open('pdb2pqr/dat/AMBER.DAT','a') as amber_file:
        amber_file.write(pdb2pqr_additional_params)
        
    
    #Handle AA.xml:
    
    #Match coordinates to the system used in pdb2pqr
    gly = Atoms(symbols = ['N','C','C','O'],
              pdbs = ['N','CA','C','O'],
              positions = np.array([
                [ 1.201, 0.847, 0.000],
                [ 0.000, 0.000, 0.000],
                [-1.250, 0.881, 0.000],
                [-2.185, 0.660,-0.784]
            ]))
    
    subset = [a.index for a in pp_nsr for i in range(4) if a.pdb == ['N','CA','C','O'][i]]
    gly.fit(pp_nsr, range(len(gly)), subset)
    
    #Rewrite AA.xml, HYDROGENS.xml and PATCHES.xml to contain nsr positions and connectivity:
    
    def append_children(element, **kwargs):
        for k,v in kwargs.iteritems():
            if not v is None:
                if not isinstance(v, list):
                    v = [v]
                for t in v:
                    child = tree.Element(k)
                    child.text = t
                    element.append(child)
                    
    def write_pretty(root, filepath):
        md = minidom.parseString(tree.tostring(root))
        output_xml = ''.join([line.strip() for line in md.toxml().splitlines()])
        md = minidom.parseString(output_xml)
        string = md.toprettyxml(indent = "   ")
        with open(filepath, "w") as xml_file:
            xml_file.write(string)
    
    aa_tree = tree.parse('pdb2pqr/dat/AA.xml')
    aa_root = aa_tree.getroot()
    patches_tree = tree.parse('pdb2pqr/dat/PATCHES.xml')
    patches_root = patches_tree.getroot()
    h_tree = tree.parse('pdb2pqr/dat/HYDROGENS.xml')
    h_root = h_tree.getroot()

    ns = pp_nsr.get_neighbours()

    for resnum, residue in pp_nsr.get_residue_dict().iteritems():

        res_el = tree.Element('residue')
        append_children(res_el, name = residue)
        
        patch_el = tree.Element('patch') 
        append_children(patch_el, name = residue + 'A', applyto = residue)
        
        class_el = tree.Element('class') 
        append_children(class_el, name = residue, opttype = "Generic")
        
        add_el = tree.Element('add')
        
        for a in pp_nsr.take(resnums = resnum, indices_in_tags = True):
            atom_el = tree.Element('atom')
            altname = a.pdb.replace("'", "*") if "'" in a.pdb else None
            append_children(atom_el, 
                            name = a.pdb, 
                            altname = altname, 
                            x = '%.3f' % a.x, 
                            y = '%.3f' % a.y, 
                            z = '%.3f' % a.z, 
                            bond = list(pp_nsr.get_pdbs()[ns[a.tag]]))
            res_el.append(atom_el)
            add_el.append(atom_el)
        aa_root.append(res_el)
        patch_el.append(add_el)
        patches_root.append(patch_el)
        h_root.append(class_el)

    write_pretty(aa_root, 'pdb2pqr/dat/AA.xml')
    write_pretty(patches_root, 'pdb2pqr/dat/PATCHES.xml')
    write_pretty(h_root, 'pdb2pqr/dat/HYDROGENS.xml')
    
    #Add hydrogens calculated from pybel to be optimised by pdb2pqr
    h_indices = pp_nsr.take('H', indices_in_tags = True).get_tags()

    #Function to find a chain of 4 heavy atoms from H
    def walk(structure,heavy_list,index):
        n = [i for i in structure.expand_selection([heavy_list[-1]], inclusive = False) if i not in heavy_list and structure[i].symbol != 'H']
        if len(n) > 0:
            try:
                heavy_list.append(n[index])
            except(IndexError):
                return False
            return True
        else:
            return False
            
    n_opt_list = []
    h_structure = Atoms()
    reduced_nsr = pp_nsr.remove('H', indices_in_tags = True).get_tags()
    for h in h_indices:
        heavy_list = pp_nsr.expand_selection([h], inclusive = False)
        index = 0
        step = 0
        while len(heavy_list) != 4:
            if not walk(pp_nsr, heavy_list, index):
                heavy_list = heavy_list[:-1]
                index += 1
            else:
                index = 0
            step += 1
            if step > 99:
                raise RuntimeError("Number of steps in heavy atom search exceeded.")

        subset = pp_nsr.take()[heavy_list+[h]]
        reduced_heavy_list = [j for i in heavy_list for j,t in enumerate(reduced_nsr) if i==t]
        protein_subset_indices = list(protein.take(residues='CHR', indices_in_tags = True)[reduced_heavy_list].get_tags())
        protein.fit(subset, protein_subset_indices, range(4))
        h_structure.append(subset[-1])
        
        #These hydrogens might clash when added to the protein so a crude dihedral optimisation is needed
        if pp_nsr[heavy_list[0]].pdb == 'N':
            n_opt_list.append(protein_subset_indices[0])
        
    protein += h_structure
    protein.calculate_neighbours()
    
    
    for nh in n_opt_list:
        hn = [h for h in protein.expand_selection([nh], inclusive=False) if protein[h].symbol == 'H'][0]
        neighbours = protein.expand_selection([nh], mode='distances', expansion = 2)
        c = [c for c in neighbours if protein[c].pdb == 'C'][0]
        ca = [ca for ca in neighbours if protein[ca].pdb == 'CA'][0]
        next_heavy = [h for h in protein.expand_selection([ca]) if h not in [nh, ca] and protein[h].symbol != 'H'][0]
        dihedral = [next_heavy, ca, nh, hn]
        r_d=[]
        for s in range(12):
            dihedral_angle = 2 * np.pi * float(s) / 12
            protein.set_dihedral(dihedral, dihedral_angle)
            distance = protein.get_distance(hn, c)
            r_d.append([distance, dihedral_angle])
        opt_dihedral=sorted(r_d)[-1][1]
        protein.set_dihedral(dihedral, opt_dihedral)
        
    
    
    os.system('cd ..')
    protein.write_pdb('protein.pdb')
    
    ph_str = '--ph-calc-method=propka --with-ph={p} '.format(p = ph) if ph else ''
    os.system('pdb2pqr/pdb2pqr -v --ff=amber --ffout=amber --nodebump {p}protein.pdb protein.pqr > pdb2pqr.log'.format(p = ph_str))
    
    if not debug:
        os.system('rm -r pdb2pqr')
        
    prot_protein = read_pqr('protein.pqr', atom_col_type = 'pdb')
    old_tags = list(reversed(protein.get_tags()))
    for a in prot_protein:
        a.tag = old_tags.pop() if a.symbol != 'H' else 0
    return prot_protein

def get_atom_parameters(atoms, filename, debug = False):
    if not filename:
        filename = atoms.info.get("name")
    if not filename:
        filename = "".join(random.choice(string.ascii_lowercase + string.digits) for i in range(6))
    dirname = filename + "_params" + os.sep
    
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        
    os.chdir(dirname)
    set_amberhome()
    atoms.write_mol2('{name}.mol2'.format(name = filename), atom_col_1 = 'pdb', atom_col_2 = 'amber')
    os.system('parmchk2 -a Y -f mol2 -i {name}.mol2 -o {name}.frcmod'.format(name = filename))
    params = get_database_parameters(os.getcwd() + os.sep + '{name}.frcmod'.format(name = filename), 'frcmod')

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
        warnings.warn('\nErrors found in %s:\n' % (filename) + warn_str)
    
    os.chdir("..")
    if not debug:
        os.system("rm -r {d}".format(d = dirname))
        
    return(params)
    
def get_database_parameters(database_name, database_type):
    
    set_amberhome()
    if database_type == 'dat':
        return get_amber_params_from_dat(database_name)
    elif database_type == 'frcmod':
        return get_amber_params_from_frcmod(database_name)
    else:
        raise RuntimeError("Database type not understood")

def set_amberhome(path = None):
    if path:
        os.environ['AMBERHOME'] = path
    else:
        if not os.environ.get('AMBERHOME'):
            if os.environ.get('CONDA_ENV_PATH'):
                try:
                    os.listdir(os.environ['CONDA_ENV_PATH'] + os.sep + "dat/leap/parm")
                    os.environ['AMBERHOME'] = os.environ['CONDA_ENV_PATH']
                except OSError:
                    raise RuntimeError("Parameter files not found in AMBERHOME")
            else:
                raise RuntimeError("AMBERHOME not set. AMBERHOME must contain parameter files in $AMBERHOME/dat/leap/parm/")

def get_oniom_model_from_indices(target,oniom_list):
    """
    Input target as an ase Atoms object
    Input indices as a list corresponding to the indices of the model region
    Returns ase Atoms object corresponding to the oniom model region
    """
    
    oniom_extraction = target.take(indices_in_tags = True)    
    indices = copy.deepcopy(oniom_list)
    
    oniom_extraction.set_calculator(Gaussian(label = 'oniom_extraction', basis = "oniom", method = "oniom(HF:UFF)=OnlyInputFiles"))
    oniom_extraction.calc.oniom_coord_params['layers'] = [indices]
    oniom_extraction.calc.oniom_coord_params['layer_mults'] = [1, 1]
    oniom_extraction.calc.set_job(time=72 ,nodes=16, memory=48000)
    oniom_extraction.calc.start(frc = True)
    
    atoms = oniom_extraction[oniom_list].take(indices_in_tags = True)
    with open("oniom_extraction.log","r") as logfile:
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
                
    
    set_link_ambers(atoms)
    return atoms
    
def set_link_ambers(atoms):
    neighbours = atoms.get_neighbours()
    for i in range(len(atoms)):
        if atoms[i].symbol == 'H' and atoms[i].amber == '':
            neighbour_set = neighbours[i]
            if len(neighbour_set) != 1: raise Exception("Wrong number of neighbours")
            neighbour = list(neighbour_set)[0]
            neighbour_type = atoms[neighbour].symbol
            if neighbour_type == "C":
                atoms[i].amber = 'HC'
            else:
                atoms[i].amber = 'H'