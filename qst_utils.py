__author__ = 'clyde'

import os

import pybel
import ase


#When running a qst3/qst2 calculation we need to match the atoms in the product/reactants/transition state,
#normally this is done through gaussview/connection but here we use babel to do the same thing algorithmically
debug=True

#algorithm adapted from http://forums.openbabel.org/match-atom-order-in-two-pdb-files-of-same-molecule-td4655723.html
def match_atom_order(lead_mol, *args):
    """sets order of atoms in the ase atoms objects present in args to the order preesnt in lead_mol.
    Assumes that all atoms objects have the same number and type of atoms"""

    ob = pybel.ob
    lead_mol.write('temp_0.xyz')
    for i, mol in enumerate(args):
        mol.write('temp_{n}.xyz'.format(n=i+1))

    omols = [pybel.readfile('xyz', 'temp_{n}.xyz'.format(n=i)).next() for i in range(len(args)+1)]

    for mol in omols:
        mol.make3D()

    query = ob.CompileMoleculeQuery(omols[0].OBMol)
    mapper = ob.OBIsomorphismMapper.GetInstance(query)
    mappings = [ob.vpairUIntUInt() for i in range(len(args))]

    for omol, mapping in zip(omols[1:], mappings):
        mapper.MapFirst(omol.OBMol, mapping)
        omol.OBMol.RenumberAtoms([x[1]+1 for x in mapping])
        print list(mapping)

    aligned_mols = []
    for i, omol in enumerate(omols):
        omol.write('xyz', 'mapped_temp_{n}.xyz'.format(n=i))
        aligned_mols.append(ase.io.read('mapped_temp_{n}.xyz'.format(n=i)))

    if not debug:
        for i in range(len(args)+1):
            os.remove('temp_{n}.xyz'.format(n=i))
            os.remove('mapped_temp_{n}_xyz'.format(n=i))

    return aligned_mols