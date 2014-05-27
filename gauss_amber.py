__author__ = 'clyde'
import os
import Biskit as B
from chem_utils import set_frag_atom_types, get_ht_ind, set_nonstd_name
from pybel import readfile
from molmod.io.xyz import XYZReader

pdb_code = '1GFL'
orig_pdb = pdb_code + '.pdb'
mod_pdb = 'mod_' + orig_pdb
prot_mod_pdb = 'prot_' + mod_pdb
prot_mod_pqr = prot_mod_pdb.replace('.pdb', '.pqr')
master_pdb = pdb_code + '_master.pdb'

#TODO algorithmically find the non_standard residues -> chimera can do this via:
#http://www.cgl.ucsf.edu/pipermail/chimera-dev/2009/000605.html
non_standard_index = [64,65,66]

gfp_raw=B.PDBModel(orig_pdb)
gfp_only = gfp_raw.takeChains([0])
nonstd = gfp_only.takeResidues(non_standard_index)
nonstd.writePdb('raw_nonstd.pdb')

#merge nonstd residues for total protein
for i in range(len(non_standard_index)-1):
    gfp_only.mergeResidues(64, name='CHR')
gfp_only.renumberResidues(start=1)
gfp_only.writePdb(mod_pdb)

clean = B.PDBCleaner( nonstd )
capped_nonstd = clean.capTerminals( nonstd, 0 )

#merge nonstd residues for nonstd part (capping code doesn't work if nonstd residues present)
for i in range(len(non_standard_index)-1):
    capped_nonstd.mergeResidues(1,name='CHR')
capped_nonstd.renumberResidues(start=64)
capped_nonstd.writePdb('raw_capped_nonstd.pdb')

#protonation that ignores the residue names work so can use babel to protonate
#WARNING Babel splits the capping residues into new residues and distributes them through out the file,
# this might be a problem in some instances - red tools doesn't seem to mind though and -out.pdb it generates puts the capping residues back where they should be
#should manually edit to make sure the caps and protons are in a sensible place/bond orders are correct)
mol = readfile("pdb", "raw_capped_nonstd.pdb").next()
mol.addh()
mol.write('pdb', 'capped_nonstd.pdb', overwrite=True)


#protonate gfp_only except for non_standard (pdb2pqr)
os.system('pdb2pqr --ff=amber --ffout=amber --nodebump {i} {o}'.format(i=mod_pdb, o= prot_mod_pqr))
os.system('babel {i} {o}'.format(i=prot_mod_pqr, o=prot_mod_pdb))

#Run Ante_red on capped_nonstd generating capped_nonstd-out.p2n
os.system('Ante_RED-1.5.pl {pdb_f}'.format(pdb_f='capped_nonstd.pdb'))

#mod .p2n to constrain charges of caps to 0
capped_nonstd = B.PDBModel('capped_nonstd-out.pdb')
ace_capping_atoms = [str(n+1) for n,r in enumerate(capped_nonstd.atoms['residue_name']) if r == 'ACE']
nme_capping_atoms = [str(n+1) for n,r in enumerate(capped_nonstd.atoms['residue_name']) if r == 'NME']

ace_remark = 'REMARK INTRA-MCC 0.0 |  ' + '  '.join(ace_capping_atoms) + ' | Remove\n'
nme_remark = 'REMARK INTRA-MCC 0.0 |  ' + '  '.join(nme_capping_atoms) + ' | Remove\n'

with open('capped_nonstd-out.p2n') as p2n_f:
    p2n_content = p2n_f.read()
split_contents = p2n_content.split('REMARK\n')

split_contents[-2] += ace_remark
split_contents[-2] += nme_remark
p2n_content = 'REMARK\n'.join(split_contents)

with open('mod_capped_nonstd-out.p2n', 'w') as p2n_f:
    p2n_f.write(p2n_content)


##use red (opt=true, qm program = gaussian) to generate mol2 file with partial charges
#run_red(mod_capped_nonstd-out.p2n)


#no nonstd residue present in file generated from pdb2pqr so now we add in the nonstd residue
#get coords of the protonated residue from the red_tools pdb output (babel mixes up the
#residues when it protonates but Ante_red fixes this)
#get amber atom names from the mol2 file we are using to make the library file
#(cannot take coords as red/gaussian end up shifting the zero coordinates)

os.system('antechamber -fi mol2 -i Mol_m1-o1-sm.mol2 -fo pdb -o mod_capped_nonstd.pdb')
amber_atom_names = B.PDBModel('mod_capped_nonstd.pdb').atoms['name']
gfp_chr = B.PDBModel('capped_nonstd-out.pdb').takeResidues([1])
gfp_chr.atoms.set('name', amber_atom_names)

prot_gfp = B.PDBModel(prot_mod_pdb)
#residues before the non_std residue
c0_res_list = range(len(prot_gfp.resList()))[0:non_standard_index[0]]
#residues after the non_std residue
c1_res_list = range(len(prot_gfp.resList()))[non_standard_index[0]:]
prot_gfp_c0 = prot_gfp.takeResidues(c0_res_list)
prot_gfp_c1 = prot_gfp.takeResidues(c1_res_list)


#build the master from the protonated standard residues and the non-standard residue (currently nonstd residue is maximally protonated - need to run propka/pdb2pqr with custom data for nonstd residue to have ph dependent protonation of nonstd residue
master_gfp = prot_gfp_c0.concat(gfp_chr).concat(prot_gfp_c1)
master_gfp.mergeChains(0)
master_gfp.mergeChains(0)
master_gfp.renumberResidues()
master_gfp.writePdb(master_pdb)
#antechamber renumbers the atoms so they are consistent (bizkit doesn't - leaves atom numbers inconsistent after the concatenation)
os.system('antechamber -fi pdb -i {m} -fo pdb -o {m}'.format(m=master_pdb))

#create set of complete list of gaff and amber atom types so we can map one to the other
os.system('antechamber -fi pdb -i {p} -fo mol2 -o {m}'.format(p=master_pdb, m=pdb_code +'_master_gaff.mol2'))
os.system('antechamber -fi pdb -i {p} -fo mol2 -o {m} -at amber'.format(p=master_pdb, m=pdb_code +'_master_amber.mol2'))
non_std_res = non_standard_index[0]
non_std_atom_indexes = [a['serial_number']-1 for a in B.PDBModel(master_pdb).resList()[non_std_res]]

#schrodinger suite utilities structconvert pdb -> mdf, computes connectivity first using standard AAs then standard heuristics
#then babel mdf -> sdf, sdf -> molmod -> neighbours
#master_sdf = SDFReader(pdb_code + '_master.sdf').next()
master_xyz = XYZReader(pdb_code + '_master.xyz').get_first_molecule()
master_xyz.set_default_graph()

#to determine the head and tail atoms of the non_standard residue we use the nearest neighbours
head_atom_ind, tail_atom_ind = get_ht_ind(master_xyz, non_std_atom_indexes)
head_atom_nm, tail_atom_nm = gfp_chr[head_atom_ind]['name'], gfp_chr[tail_atom_ind]['name']

#antechamber mol2 -> mol2 to generate a mol2 file that uses amber names fails because for the backbone N and C
#because we have cut bonds and antechamber assumes that means they are double bonded!
#so we need to take the names from the full protein
os.system('antechamber -fi mol2 -i Mol_m1-o1-sm.mol2 -fo mol2 -o amb_Mol_m1-o1-sm.mol2')
set_frag_atom_types('amb_Mol_m1-o1-sm.mol2', 'mod_capped_nonstd.mol2', master_pdb, atom_types='amber')
set_nonstd_name('mod_capped_nonstd.mol2', 'CHR')
#use parmchk2 on mol2 to generate .frcmod file
os.system('parmchk2 -a Y -f mol2 -i mod_capped_nonstd.mol2 -o mod_capped_nonstd.frcmod')

#if we use gaff atoms this defines a gaff forcefield for the non standard residue and specifies gaff atom
#types along with the parameters defining their interaction however the rest of the protein is defined by
#an amber forcefield so we need to define parameters that define the interaction between the terminal n and
#c gaff atoms and the amber atom types that it can join too. To do this we can use the standard amber parameters.

#gaff_atoms, gaff_bonds = mol2_parse(pdb_code +'_master_gaff.mol2')
#amber_atoms, amber_bonds = mol2_parse(pdb_code +'_master_amber.mol2')
#get bonds,angles and dihedrals required for atoms spanning the nonstandard residues and the standard residues
#master_lookup = get_bad(master_xyz, gaff_atoms, amber_atoms, non_std_atom_indexes)
#bonds, angles, dihedrals = master_lookup['mixed']
#amber_bonds, amber_angles, amber_dihedrals = master_lookup['amber']
#extra_params = get_amber_params([amber_bonds, amber_angles, amber_dihedrals)
#add_params('mod_capped_nonstd.frcmod', extra_params)

"""Could not find bond parameter for: C - n
Could not find bond parameter for: c - N
Building angle parameters.
Could not find angle parameter: O - C - n
Could not find angle parameter: C - n - hn
Could not find angle parameter: C - n - c3
Could not find angle parameter: CT - C - n
Could not find angle parameter: o - c - N
Could not find angle parameter: c - N - H
Could not find angle parameter: c - N - CT
Could not find angle parameter: c3 - c - N
Building proper torsion parameters.
 ** No torsion terms for  O-C-n-hn
 ** No torsion terms for  O-C-n-c3
 ** No torsion terms for  CT-C-n-hn
 ** No torsion terms for  CT-C-n-c3
 ** No torsion terms for  o-c-N-H
 ** No torsion terms for  o-c-N-CT
 ** No torsion terms for  c3-c-N-H
 ** No torsion terms for  c3-c-N-CT"""



#generate non_std.in and master.in
#run tleap using non_std.in to generate library, file

#set CHR.1 connect0 CHR.1.{head_atom}
#set CHR.1 connect1 CHR.1.{tail_atom}
#set CHR.1 restype protein
#set CHR.1 name "CHR"

tleap_1_content = """source leaprc.ff99SB
source leaprc.gaff
loadamberparams mod_capped_nonstd.frcmod
CHR = loadmol2 mod_capped_nonstd.mol2
set CHR head CHR.1.{head_atom}
set CHR tail CHR.1.{tail_atom}
set CHR name "CHR"
check CHR
saveoff CHR CHR.lib
saveamberparm CHR chr.prmtop chr.inpcrd
quit """.format(head_atom=head_atom_nm, tail_atom=tail_atom_nm)

with open('non_std.in', 'w') as tleap_in_f:
    tleap_in_f.write(tleap_1_content)

os.system('tleap -f {inp_fn}'.format(inp_fn='non_std.in'))

#run tleap using master.in to generate final params for gfp_master
tleap_2_content = """source leaprc.ff99SB
source leaprc.gaff
loadamberparams mod_capped_nonstd.frcmod
loadoff CHR.lib
complex = loadpdb {mast}
check complex
saveamberparm complex 1GFL_chr.prmtop 1GFL_chr.inpcrd
savepdb complex 1GFL_chr.pdb
quit """.format(mast=master_pdb)

with open(pdb_code + '.in', 'w') as tleap_in_f:
    tleap_in_f.write(tleap_2_content)

os.system('tleap -f {inp_fn}'.format(inp_fn=pdb_code + '.in'))