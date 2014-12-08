__author__ = 'clyde'

from ase import Atoms
from molmod.io import FCHKFile

#Gaussian constant determining no. angstrom per bohr
#see http://www.gaussian.com/g_tech/g_ur/k_constants.htm
ang_per_bohr = 0.5291772086


class FCHK(FCHKFile):
    @property
    def irc_no_frames(self):
        return self.fields.get('IRC Number of geometries', [0])[0]

    @property
    def irc_energies(self):
        bare_energy_steps = self.fields['IRC point       1 Results for each geome']
        frame_energies = [bare_energy_steps[i] for i in range(0, len(bare_energy_steps), 2)]
        return frame_energies

    @property
    def irc_steps(self):
        bare_energy_steps = self.fields['IRC point       1 Results for each geome']
        frame_steps = [bare_energy_steps[i] for i in range(1, len(bare_energy_steps), 2)]
        return frame_steps

    @property
    def atomic_numbers(self):
        atomic_numbers = self.fields['Atomic numbers']
        return atomic_numbers

    @property
    def irc_frames(self):
        no_atom_coords = self.fields['IRC Num geometry variables']

        bare_coords= self.fields['IRC point       1 Geometries']
        bare_frame_coords = [bare_coords[f*no_atom_coords:(f+1)*no_atom_coords] for f in range(self.irc_no_frames)]
        frame_coords = [[c[a*3:(a+1)*3]*ang_per_bohr for a in range(no_atom_coords/3)] for c in bare_frame_coords]

        irc_frames = [Atoms(numbers=self.atomic_numbers, positions=coords) for coords in frame_coords]
        return irc_frames

    @property
    def final_geom(self):
        return self.molecule.coordinates

    def sort_frames(self):
        #mid_index=self.irc_no_frames/2
        #first_frames, second_frames = self.irc_frames[mid_index+1:], self.irc_frames[0:mid_index+1]
        #first_frame_es, second_frame_es = self.irc_energies[mid_index+1:], self.irc_energies[0:mid_index+1]
        #first_frame_steps, second_frame_steps = self.irc_steps[mid_index+1:], self.irc_steps[0:mid_index+1]

        #first_frames.reverse(), first_frame_es.reverse(), first_frame_steps.reverse()

        #irc_frames = first_frames + second_frames
        #frame_energies = first_frame_es + second_frame_es
        #frame_steps = first_frame_steps + second_frame_steps
        return