__author__ = 'clyde'

import warnings

import numpy as np

from my_cc_utils.data_utils import unwind


#utilities
normalize = lambda e: e/np.linalg.norm(e)


class Vibration(object):
    def __init__(self, freq, disp):
        self._freq = freq
        self._disp = unwind(disp)

    @property
    def freq(self):
        return self._freq

    @property
    def disp(self):
        return self._disp

    def overlap(self, other):
        "Returns a normalized absolute overlap between two modes based on their displacements"
        norm_disp = normalize(self.disp)
        norm_other_disp = normalize(other.disp)

        return abs(np.dot(norm_disp, norm_other_disp))


class Vibrations(object):
    def __init__(self, list_vibrations):
        self._vibs = list_vibrations

    def __getitem__(self, key):
        return self.vibs[key]

    def add(self,vib):
        self.vibs.append(vib)

    def __iter__(self):
        for vib in self.vibs:
            yield vib

    def __len__(self):
        return len(self.vibs)

    @property
    def vibs(self):
        return self._vibs

    @property
    def freqs(self):
        return [vib.freq for vib in self.vibs]

    @property
    def disps(self):
        return [vib.disp for vib in self.vibs]

    def sort_by_index(self, sort_indices):
        """Sorts the vibrations according to a list of indices"""

        temp_enum_vibs = zip(self.vibs, sort_indices)
        temp_enum_vibs.sort(key=lambda e: e[1])
        self.vibs, sorted_inds = zip(*temp_enum_vibs)


    def alignment(self, vec):
        """Returns a a list of overlaps (defined between 0 and 1) for each vibration amplitude and the vector given"""
        als = []
        norm_vec = normalize(vec)

        for i in range(len(self.vibs)):
            norm_mode_disp = normalize(self.vibs[i].disp)
            als.append(abs(np.dot(norm_mode_disp, norm_vec)))

        return als

    def alignments(self, list_vecs):
        """Returns a list of overlaps (defined between 0 and 1) for each vibration amplitudes
        and each of the list of vectors given"""

        if len(self) != len(list_vecs):
            raise(RuntimeError('Alignment requires two objects of equal size'))

        als = []
        for i in range(len(self.vibs)):
            norm_mode_disp = normalize(self.vibs[i].disp)
            norm_vec = normalize(list_vecs[i])
            als.append(abs(np.dot(norm_mode_disp, norm_vec)))

        return als

    def mode_match(self, other):
        """Returns a list of tuples each of which consists of two indices and a score.

        The indices map the modes in this object to the modes in the other object
        whilst the score represents the overlap between the identified modes.

        I.e. if self = [m1,m2,m3,m4] and other = [m2,m1,m3,m4],
        self.mode_match(other) returns [(0,1),(1,0),(2,2),(3,3)]"""

        if len(self) != len(other):
            raise(RuntimeError('Cannot track different numbers of modes'))

        overlap_data = []
        for self_ind,mode in enumerate(self):
            overlaps = [mode.overlap(other_mode) for other_mode in other]
            sorted_overlaps = sorted(overlaps)

            if sorted_overlaps[-1] == sorted_overlaps[-2]:
                warnings.warn("Mutiple equivalent best matches for mode={i} freq={f}".format(i=self_ind, f=mode.freq))

            max_overlap = sorted_overlaps[-1]
            matching_other_ind = overlaps.index(max_overlap)

            overlap_data.append((self_ind,matching_other_ind, max_overlap))

        low_overlaps = [(s_ind,o_ind,overlap) for (s_ind,o_ind, overlap) in overlap_data if overlap<0.8]

        if any(low_overlaps):
            warnings.warn('Overlaps less than 80% detected: {l_o}'.format(l_o=low_overlaps))

        return overlap_data

def get_modes(ase_frame):
    """Returns a Vibrations object contains the modes present in a gaussian calculation associated with an ASE object"""
    freqs = ase_frame.calc.max_data['vibfreqs']
    disps = ase_frame.calc.max_data['vibdisps']
    frame_modes = Vibrations([Vibration(*vib) for vib in zip(freqs,disps)])
    return frame_modes

def sort_by_mode_id(list_vibs):
    """Takes a list of Vibrations objects and aligns them so that modes at index[i] are maximally overlapped.

    I.e. [m1,m2,m3,m4], [m1',m2',m4',m3'], [m1'',m2'',m3'',m4''] becomes
         [m1,m2,m3,m4], [m1',m2',m3',m4'], [m1'',m2'',m3'',m4'']

    Modes are matches from the left to the right, i.e. the first Vibrations object sets the order.

    NB we are assuming that each Vibrations object only differs from it's neighbouring objkects in the list by a
    small perturbation. Inconsistent behaviour will be seen if we attempt to sort very different sets of modes"""

    sorted_list_vibs = [list_vibs[0]]
    for i, modes in enumerate(list_vibs):
        if i != 0:
            map_ind = zip(*modes.mode_match(previous_modes))[1]
            sorted_modes = Vibrations([modes[i] for i in map_ind])
            sorted_list_vibs.append(sorted_modes)
            previous_modes = sorted_modes
        else:
            previous_modes = modes
    return sorted_list_vibs