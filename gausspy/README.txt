

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
|        +<-+SR{prot} |     |NSR{prot}+-----+
|        |  +----+----+     +----+----+     |
|        |       |               |          |
|        v       |               |          v
|     Extract    +------>+<------+  +------>+
|       MR               |          |       |
|        |               |          |  +====v====+
|   +----v----+    +-----v-----+    |  | Params  |
|   | MR{prot}|    |Chain{prot}+----+  |{SR, NSR}|
|   +----+----+    +-----+-----+       +====+====+
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
* SR:    Standard residues     *
* NSR:   Non-standard residues *
*                              *
********************************

1.2 Overview

protonate() takes in a Chain as an ASE Protein Atoms object along with the
NSR indices (residue numbers) and separately protonates the NSRs and SRs.
It returns the protonated chain.

1.3 Usage

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

1.4 Example

The following code can be used to protonate and perform an opt calculation on
 1GFL.pdb


from ase_extensions.protein_utils import *
from gausspy.oniom_amber_utils import *
from gausspy import Gaussian

protein=read_pdb('1GFL.pdb')
chain_a=protein.take(chains='A',indices_in_tags=True).remove(residues='HOH')

#NSR is in residues 65-67
non_standard_indices=[65,66,67]
#These are the atom numbers corresponding to the model region in the initial
#chain. The tags will be remembered during protonation so extraction is trivial
old_oniom_list=[479]+range(482,495)
master_atoms=protonate(atoms=chain_a,
                        non_standard_indices=non_standard_indices,
                        use_red=False,
                        ph=7,
                        force_capping=True,
                        non_standard_residue_name='CHR',
                        debug_level=4)

nsr=master_atoms.take(residues='CHR')
#The new indices of the model region are calculated here
oniom_list=list(master_atoms.take(tags=old_oniom_list,
                                       indices_in_tags=True).get_tags())
#Here we find the neighbouring hydrogens to add them to the model region
neighbour_list=master_atoms.expand_selection(oniom_list,
                                             mode='bonds',
                                             expansion=1,
                                             inclusive=False)
h_list=[n for n in neighbour_list if master_atoms[n].symbol=='H']
oniom_list+=h_list
#This function returns an atoms object corresponding to the model region and 
#the link atoms.
model_region=get_oniom_model_from_indices(master_atoms,sorted(oniom_list))

#Parameters are calculated here
parm99_params=get_database_parameters('parm99','dat')
ff99SB_params=get_database_parameters('ff99SB','frcmod')
nsr_params=get_atom_parameters(nsr,'nsr',True)
model_params=get_atom_parameters(model_region,'model',True)

#Any additional parameters can be calculated this way, in case the Amber 
#calculation fails.
additional_params=get_atom_parameters(master_atoms.take(resnums=range(68,71)),
                                                        'additional')

#An initial Amber calculation (optional) is used to test the parameters and
#optimise the structure
master_atoms.set_calculator(Gaussian(label="amber",
                                     opt="opt",
                                     method="amber=softfirst",
                                     geom="connectivity"))
master_atoms.calc.set_job(time=72,nodes=4,memory=4000)
master_atoms.calc.set_amber_params([model_params,
                                    nsr_params,
                                    parm99_params,
                                    ff99SB_params,
                                    additional_params])

master_atoms.calc.start(frc=True)

#Here, we can copy the positions of the atoms from the previous calculations
master_atoms.positions = copy.deepcopy(master_atoms.calc.read_atoms().positions)

#The ONIOM calculation is performed here
master_atoms.set_calculator(Gaussian(label="test",
                                     basis="oniom",
                                     opt="opt",
                                     method="oniom(B3LYP/6-31G(d):amber=softfirst)",
                                     geom="connectivity"))
master_atoms.calc.oniom_coord_params['layers']=[oniom_list]
master_atoms.calc.set_job(time=72,nodes=4,memory=4000)
master_atoms.calc.set_amber_params([model_params,additional_params,nsr_params,parm99_params,ff99SB_params])

master_atoms.calc.start(frc=True)
