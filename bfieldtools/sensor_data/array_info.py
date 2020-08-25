# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 14:22:27 2020

@author: Antti
"""


def extract_array_info(sensor_type="mag"):
    """
    Extract array info from mne
    
    For this mne with sample_data needs to be installed
    
    sensor_type value must be one of 
    ['grad', 'mag', 'planar1', 'planar2'] or True (all channels)

    Returns
    -------
    dictionary of transformation matrices

    """
    from mne import read_evokeds
    from mne.datasets import sample

    data_path = sample.data_path()
    fname = data_path + "/MEG/sample/sample_audvis-ave.fif"
    # Reading
    condition = "Left Auditory"
    evoked = read_evokeds(fname, condition=condition, verbose=False)
    evoked.pick_types(meg=sensor_type)

    def loc2mat(loc):
        mat = np.eye(4)
        mat[:3, 3] = loc[:3]
        mat[:3, 0] = loc[3:6]
        mat[:3, 1] = loc[6:9]
        mat[:3, 2] = loc[9:]

        return mat

    names = [ch["ch_name"].replace(" ", "") for ch in evoked.info["chs"]]
    locs = [ch["loc"] for ch in evoked.info["chs"]]
    mats = [loc2mat(loc) for loc in locs]

    return dict(zip(names, mats))
