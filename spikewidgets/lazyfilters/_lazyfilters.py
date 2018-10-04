from abc import ABC, abstractmethod
import spikeinterface as si
import numpy as np
from scipy import special

class LazyFilterRecording(si.InputExtractor):
    def __init__(self, *, recording,chunk_size=10000):
        si.InputExtractor.__init__(self)
        self._recording=recording
        self._chunk_size=chunk_size
        self._filtered_chunks=dict()
        self.getNumChannels=recording.getNumChannels
        
    def getNumChannels(self):
        return self._recording.getNumChannels()
    
    def getNumFrames(self):
        return self._recording.getNumFrames()
    
    def getSamplingFrequency(self):
        return self._recording.getSamplingFrequency()
        
    def getRawTraces(self, start_frame=None, end_frame=None, channel_ids=None):
        if start_frame is None:
            start_frame=0
        if end_frame is None:
            end_frame=self.getNumFrames()
        if channel_ids is None:
            channel_ids=list(range(self.getNumChannels()))
        ich1=int(start_frame/self._chunk_size)
        ich2=int((end_frame-1)/self._chunk_size)
        filtered_chunk_list=[]
        for ich in range(ich1,ich2+1):
            filtered_chunk0=self._get_filtered_chunk(ich)
            if ich==ich1:
                start0=start_frame-ich*self._chunk_size
            else:
                start0=0
            if ich==ich2:
                end0=end_frame-ich*self._chunk_size
            else:
                end0=self._chunk_size
            filtered_chunk_list.append(filtered_chunk0[channel_ids,start0:end0])
        return np.concatenate(filtered_chunk_list,axis=1)
    
    @abstractmethod
    def filterChunk(self,*,start_frame,end_frame):
        raise NotImplementedError('filterChunk not implemented')
    
    def _get_filtered_chunk(self, ind):
        code=str(ind)
        if code not in self._filtered_chunks:
            start0=ind*self._chunk_size
            end0=(ind+1)*self._chunk_size
            self._filtered_chunks[code]=self.filterChunk(start_frame=start0,end_frame=end0)
        return self._filtered_chunks[code]
    
    def getChannelInfo(self, channel_id):
        return self._recording.getChannelInfo(channel_id=channel_id)
    
class BandpassFilterRecording(LazyFilterRecording):
    def __init__(self, *, recording, freq_min, freq_max, freq_wid):
        LazyFilterRecording.__init__(self, recording=recording, chunk_size=3000*10)
        self._recording=recording
        self._freq_min=freq_min
        self._freq_max=freq_max
        self._freq_wid=freq_wid
    def filterChunk(self,*,start_frame,end_frame):
        padding=3000
        i1=start_frame-padding
        i2=end_frame+padding
        padded_chunk=self._read_chunk(i1,i2)
        filtered_padded_chunk=self._do_filter(padded_chunk)
        return filtered_padded_chunk[:,start_frame-i1:end_frame-i1]
    def _create_filter_kernel(self,N,samplerate,freq_min,freq_max,freq_wid=1000):
        # Matches ahb's code /matlab/processors/ms_bandpass_filter.m
        # improved ahb, changing tanh to erf, correct -3dB pts  6/14/16    
        T = N / samplerate # total time
        df = 1 / T # frequency grid
        relwid = 3.0; # relative bottom-end roll-off width param, kills low freqs by factor 1e-5.

        k_inds=np.arange(0,N)
        k_inds=np.where(k_inds<=(N+1)/2,k_inds,k_inds-N)

        fgrid=df*k_inds
        absf=np.abs(fgrid)

        val=np.ones(fgrid.shape)
        if freq_min!=0:
            val=val*(1 + special.erf(relwid * (absf - freq_min) / freq_min)) / 2
            val=np.where(np.abs(k_inds)<0.1,0,val) # kill DC part exactly
        if freq_max!=0:
            val=val*(1 - special.erf((absf - freq_max) / freq_wid)) / 2;
        val=np.sqrt(val) # note sqrt of filter func to apply to spectral intensity not ampl
        return val
    def _do_filter(self,chunk):
        samplerate=self._recording.getSamplingFrequency()
        M=chunk.shape[0]
        chunk2=chunk
        # Subtract off the mean of each channel unless we are doing only a low-pass filter
        #if self._freq_min!=0:
        #    for m in range(M):
        #        chunk2[m,:]=chunk2[m,:]-np.mean(chunk2[m,:])
        # Do the actual filtering with a DFT with real input
        chunk_fft=np.fft.rfft(chunk2) 
        kernel=self._create_filter_kernel(
            chunk2.shape[1],
            samplerate,
            self._freq_min,self._freq_max,self._freq_wid
        )
        kernel=kernel[0:chunk_fft.shape[1]] # because this is the DFT of real data
        chunk_fft=chunk_fft*np.tile(kernel,(M,1))
        chunk_filtered=np.fft.irfft(chunk_fft)
        return chunk_filtered
    def _read_chunk(self,i1,i2):
        M=self._recording.getNumChannels()
        N=self._recording.getNumFrames()
        if i1<0:
            i1b=0
        else:
            i1b=i1
        if i2>N:
            i2b=N
        else:
            i2b=i2
        ret=np.zeros((M,i2-i1))
        ret[:,i1b-i1:i2b-i1]=self._recording.getRawTraces(start_frame=i1b,end_frame=i2b)
        return ret

def bandpass_filter(recording,freq_min=300, freq_max=6000, freq_wid=1000):
    return BandpassFilterRecording(
        recording=recording,
        freq_min=freq_min,
        freq_max=freq_max,
        freq_wid=freq_wid
    )
    
    