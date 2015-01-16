# Copyright (C) 2010 by CAMd, DTU
# Please see the accompanying LICENSE file for further information.

# This file is taken (almost) verbatim from CMR with D. Landis agreement

# modified 12/12/2012 by C.Fare to expand the reader by including cclib's gaussian parser
import cclib, logging

FIELD_SEPARATOR="\\"
PARA_START="\n\n"
PARA_END="\\\\@"

names = ["", "", "Computer_system", "Type_of_run", "Method", "Basis_set", "Chemical_formula", "Person", "Date", "", "", "Root", "", "Title", ""] #[Charge,Multi]
#        0   1            2               3           4            5              6               7        8    9  10  11  12    13     14
names_compact = ["", "", "Computer_system", "Type_of_run", "Method", "Basis_set", "Chemical_formula", "Person", "Date", "", "", "", "", "Title", ""] #[Charge,Multi]
#                0   1            2               3           4            5              6               7        8    9  10  11  12    13     14

charge_multiplicity = 15


class GaussianReader:

    def auto_type(self, data):
        """ tries to determine type"""
        try:
            return float(data)
        except ValueError:
            pass

        try:
            ds = data.split(",")
            array = []

            for d in ds:
                array.append(float(d))

            return array
        except ValueError:
            pass

        return data


    def __init__(self, filename, inc_cclib=False):
        """filename is optional; if not set, use parse to set the content"""
        if not filename is None:
            fin = file(filename)
            content = fin.read()
            fin.close()
            #handles the case that users used windows after the calculation:
            content = content.replace("\r\n", "\n")
            #handles the case that users used a development version of Gaussian
            content=content.replace("\n The archive entry for this job was punched.","")

            if inc_cclib:
                self.cclib_obj = cclib.parser.ccopen(filename, loglevel=logging.CRITICAL).parse()
                self.cclib_obj.listify()
                self.cclib_data = self.cclib_obj.getattributes()
            self.parse(content)

    def parse(self, content):
        #the output for atomic coords depends on whether we are doing an optimisation (single point calcs take the form 'atom, flag, x, y ,z' rather than simply 'atom, x, y, z'
        opt = 'opt' in content

        from ase.data import atomic_numbers
        self.data = []
        temp_items = content.split(PARA_START)
        seq_count = 0
        for i in temp_items:
            i=i.replace("\n ", "")
            if i.endswith(PARA_END):
                i = i.replace(PARA_END, "")
                i = i.split(FIELD_SEPARATOR)

                new_dict = {}
                self.data.append(new_dict)

                new_dict["Sequence number"] = seq_count
                seq_count += 1
                for pos in range(len(names)):
                    if names[pos]!="":
                        new_dict[names[pos]] = self.auto_type(i[pos])

                chm = i[charge_multiplicity].split(",")
                new_dict["Charge"]       = int(chm[0])
                new_dict["Multiplicity"] = int(chm[1])

                #Read atoms
                atoms = []
                positions = []
                position = charge_multiplicity+1

                if not opt:
                    atom_flags = []
                    while position<len(i) and i[position]!="":
                        s = i[position].split(",")
                        atoms.append(atomic_numbers.get(s[0], 0))

                        positions.append([float(s[-3]), float(s[-2]), float(s[-1])])
                        position = position + 1
                        atom_flags.append(int(s[1]))

                    new_dict["Atom_flags"]=atom_flags
                    new_dict["Atomic_numbers"]=atoms
                    new_dict["Positions"]=positions

                else:
                    while position<len(i) and i[position]!="":
                        s = i[position].split(",")
                        atoms.append(atomic_numbers.get(s[0], 0))

                        positions.append([float(s[-3]), float(s[-2]), float(s[-1])])
                        position = position + 1

                    new_dict["Atomic_numbers"]=atoms
                    new_dict["Positions"]=positions

                #read orbital swaps
                orb_swaps = []
                if 'alter' in new_dict.get('Root',""):
                    position +=1
                    while position<len(i) and i[position]!="":
                        s = i[position].split(",")
                        orb_swaps.append(s)
                        position+=1
                    new_dict["Orbital_swaps"]=orb_swaps

                #Read more variables
                position +=1
                while position<len(i) and i[position]!="":
                    s = i[position].split("=")
                    if len(s)==2:
                        new_dict[s[0]]=self.auto_type(s[1])
                    else:
                        print "Warning: unexpected input ",s
                    position = position + 1

    def __iter__(self):
        """returns an iterator that iterates over all keywords"""
        return self.data.__iter__()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, pos):
        return self.data[pos]
