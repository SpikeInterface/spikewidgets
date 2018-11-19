import spikeinterface as si
import numpy as np


def real(name='franklab_tetrode', download=True):
    if name == 'franklab_tetrode':
        dsdir = 'kbucket://b5ecdf1474c5/datasets/neuron_paper/franklab_tetrode'
        IX = si.MdaRecordingExtractor(dataset_directory=dsdir, download=download)
        return (IX, None)
    else:
        raise Exception('Unrecognized name for real dataset: ' + name)
