

**************************************************
********************PROTONATION*******************
**************************************************

1.1 FLOW CHART

     +-------+      +---------+      +-------+
     |  MR   |      | Initial |      |  NSR  |
     |Indices|      |PDB  File|      |Indices|
     +---+---+      +----+----+      +---+---+
         |               |               |
         |               v               |
         |         Extract Chain         |
         |          Remove  H2O          |
         |               |               |
         |            +--v--+            |
         |            |Chain|            |
         |            +--+--+            |
         |               v               |
         |               +<--------------+
         |               |
         |               v
         |            Extract
         |              NSR
         |               |
         |       +------<+>------+
         |       |               |
         |     +-v-+           +-v-+
         |     |SR |           |NSR|
         |     +-+-+           +-+-+
         |       |               |
         |       v               v
         |   Protonate       Protonate
         | with  pdb2pqr    with  Pybel
         |       |               |
+-------<+       v               v
|        v  +----+----+     +----+----+
|        |  |SR{prot} |     |NSR{prot}+-----+
|        |  +----+----+     +----+----+     |
|        |       |               |          |
|        |       |               |          |
|        |       +------>+<------+          |
|        |               |                  |
|        |               |                  |
|        v         +-----v-----+            v
|     Extract<-----+Chain{prot}+----------->+
|       MR         +-----+-----+            |
|        |               |             +====v====+
|   +----v----+          |             | Params  |
|   | MR{prot}|          |             |{SR, NSR}|
|   +----+----+          |             +====+====+
|        |               v                  |
|        |               +<-----------------+
|        |               |
|   +====v====+    +*****v*****+
|   | Params  |    *   AMBER   *
|   |  {MR}   |    *Calculation*
|   +====+====+    +*****+*****+
|        v               v
+------->+-------------->+
                         |
                   +*****v*****+
                   *   ONIOM   *
                   *Calculation*
                   +***********+


********************************
*          = KEY =             *
*                              *
* MR:    Model Region          *
* SR:    Standard Residues     *
* NSR:   Non-Standard Residues *
*                              *
********************************

1.2 OVERVIEW

protonate() takes in a Chain as an ASE Protein Atoms object along with the
NSR indices (residue numbers) and separately protonates the NSRs and SRs.
It returns the protonated chain.

1.3 USAGE

atoms:
 The chain to be protonated, consisting of standard and non-standard residues.

non_standard_indices:
 The residue numbers of the NSRs.

non_standard_residue_name:
 The name to rename all the NSRs when they are merged. 
 Recommended to be 3 characters. The default is 'CHR'.

nsr_prot_function:
 The protonation function to protonate the NSR. This function must be in the
  form:
 function('nsr_capped.pdb', debug)
 where nsr_capped.pdb is the filename corresponding to a capped, unprotonated 
  NSR.
 Debug will be set to true if protonate() has debug_level set to >2.
 The function must return the path of the capped, protonated NSR.
 It can be a dummy function that simply returns the path of an already
  protonated NSR.

use_red:
 This specifies whether to run RED III to calculate partial charges on the NSR.
 The RESP calculations is the most time-consuming part of protonate(), and so
 if this calculation was successful, set this to false to prevent rerunning it.

ph:
 The ph to use for the SR protonation.

force_capping:
 Loosens the criteria for the search for locations to cap on the NSR.
 This should be used if PDB names for the backbone in the NSR include numbers.

debug_level:
 Used to suppress warnings if set to 0, or add information to stdout.
 More detail in the function's documentation.

logfile:
 The location of the file to send logging information to.

pdb2pqr_path:
 Where the pdb2pqr program is installed.

1.4 EXAMPLE

The following code can be used to protonate and perform an opt calculation on
 1GFL.pdb


from ase_extensions.protein_utils import *
from gausspy.oniom_amber_utils import *
from gausspy import Gaussian

protein=read_pdb('1GFL.pdb')
chain_a=protein.take(chains='A',indices_in_tags=True).remove(residues='HOH')

#NSR is in residues 65-67
#In 1GFL.pdb they are given standard residue names which will break pdb2pqr,
# so their names are changed
for a in chain_a:
    if a.resnum in [65,66,67]:
        a.residue = 'CHR'

p = Protonation_Parameterisation(chain_a)
#These are the atom numbers corresponding to the model region in the initial
#chain. The tags will be remembered during protonation so extraction is trivial
p.params["Initial Model Atom Numbers"] = [479] + range(482,495)
p.auto_protonate()
p.get_nsr_charges()
p.get_model_region()
p.auto_parameterise()
model_list = p.get_final_model_atomnos()

#An initial Amber calculation (optional) is used to test the parameters.
master_atoms.set_calculator(Gaussian(label="amber",
                                     method="amber=softfirst",
                                     geom="connectivity"))
master_atoms.calc.set_job(time=72,nodes=4,memory=4000)
master_atoms.calc.set_amber_params(p.mm_parameters.values())
master_atoms.calc.start(frc=True)

#The ONIOM calculation is performed here
master_atoms.set_calculator(Gaussian(label="test",
                                     basis="oniom",
                                     opt="opt",
                                     method="oniom(B3LYP/6-31G(d):amber=softfirst)",
                                     geom="connectivity"))
master_atoms.calc.oniom_coord_params['layers']=[model_list]
master_atoms.calc.set_job(time=72,nodes=4,memory=4000)
master_atoms.calc.set_amber_params(p.mm_parameters.values())
master_atoms.calc.start(frc=True)

1.5 CHROMOPHORE DETECTION

A method for finding a candidate for the chromophore has been included. It searches 
 the entire protein for sp2 and aromatic atoms and groups them by their neighbouring
 atoms. The largest groups are proposed as candidates.

It can be used by calling 

p.find_chromophore()

before calling

p.get_model_region()

Multiple chromophores can be included in the calculation by supplying a list of indices
 to the index parameter. 

1.6 MODIFYING CHROMOPHORE CHARGE

Chromophore charges can be determined as followed:

sum(p.model_region.get_amber_charges())

Under normal conditions this is set to 0. There are two ways of changing this:

1: Modify the protonation function. This can even be a dummy function that returns a 
   previously constucted chromophore. The requirement is that it must be a Protein 
   Atoms Object with amber and pdb types calculated.

2: Modify the protonated_nsr that contains the chromophore. Run 
   p.auto_protonate() to extract the protonated NSR first. Then perform
   necessary modifiations to p.protonated_nsrs and p.protonated atoms. 
   For example:

del p.protonated_nsrs["CRA"][36]
del p.protonated_atoms[988]

   deletes atom 36 inside NSR 'CRA' (here this is an alcoholic proton)

p._hidden_params["NSR Charges Multiplicities"]["CRA"] = [-1, 1]

   sets the total charge to -1 in the RESP calculation.
