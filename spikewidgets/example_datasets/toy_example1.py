import spikeinterface as si
import numpy as np
from .synthesize_random_waveforms import synthesize_random_waveforms
from .synthesize_random_firings import synthesize_random_firings
from .synthesize_timeseries import synthesize_timeseries
import numpy as np

def toy_example1():
    num_channels=4
    duration=10
    samplerate=30000
    K=10
    upsamplefac=13

    waveforms,geom=synthesize_random_waveforms(K=K,M=num_channels,average_peak_amplitude=-100,upsamplefac=upsamplefac)
    times,labels=synthesize_random_firings(K=K,duration=duration,samplerate=samplerate)
    labels=labels.astype(np.int64)
    OX=si.NumpyOutputExtractor()
    OX.setTimesLabels(times,labels)
    X=synthesize_timeseries(output_extractor=OX,waveforms=waveforms,noise_level=10,samplerate=samplerate,duration=duration,waveform_upsamplefac=upsamplefac)

    IX=si.NumpyInputExtractor(timeseries=X,samplerate=samplerate,geom=geom)
    return (IX,OX)