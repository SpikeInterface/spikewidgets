from scipy.optimize import linear_sum_assignment
import spikeinterface as si
import numpy as np

def match_sortings(sorting,reference_sorting,delta_tp=20):
    unit_ids_ref=reference_sorting.getUnitIds()
    unit_ids=sorting.getUnitIds()
    scores=np.zeros((len(unit_ids_ref),len(unit_ids)))
    for i_ref,unit_ref in enumerate(unit_ids_ref):
        times_ref=reference_sorting.getUnitSpikeTrain(unit_ref)
        for i, unit in enumerate(unit_ids):
            times=sorting.getUnitSpikeTrain(unit)
            num_matching_events=count_matching_events(times,times_ref,delta=delta_tp)
            score=num_matching_events/(len(times)+len(times_ref)-num_matching_events)
            scores[i_ref,i]=score
    [ref_inds,inds]=linear_sum_assignment(-scores)
    ref_inds=list(ref_inds)
    inds=list(inds)
    unit_map=dict()
    k=np.max(unit_ids_ref)+1
    for i, unit in enumerate(unit_ids):
        if i in inds:
            aa=inds.index(i)
            unit_map[unit]=unit_ids_ref[aa]
        else:
            unit_map[unit]=k
            k=k+1
    return MappedSortingExtractor(sorting,unit_map)
    
class MappedSortingExtractor(si.SortingExtractor):
    def __init__(self,sorting,unit_map):
        si.SortingExtractor.__init__(self)
        self._sorting=sorting
        self._unit_map=unit_map
        self._reverse_map=dict()
        for key in unit_map:
            self._reverse_map[unit_map[key]]=key
        units=sorting.getUnitIds()
        self._unit_ids=list(np.sort([self._unit_map[unit] for unit in units]))
        
    def getUnitIds(self):
        return self._unit_ids
    
    def getUnitSpikeTrain(self,unit_id):
        unit2=self._reverse_map[unit_id]
        return self._sorting.getUnitSpikeTrain(unit2)
    
def count_matching_events(times1,times2,delta=20):
    times_concat=np.concatenate((times1,times2))
    membership=np.concatenate((np.ones(times1.shape)*1,np.ones(times2.shape)*2))
    indices=times_concat.argsort()
    times_concat_sorted=times_concat[indices]
    membership_sorted=membership[indices]
    diffs=times_concat_sorted[1:]-times_concat_sorted[:-1]
    inds=np.where((diffs<=delta)&(membership_sorted[0:-1]!=membership_sorted[1:]))[0]
    if (len(inds)==0):
        return 0
    inds2=np.where(inds[:-1]+1!=inds[1:])[0]
    return len(inds2)+1