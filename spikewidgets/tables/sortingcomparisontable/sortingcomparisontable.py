import pandas as pd
from IPython.display import HTML

class SortingComparisonTable():
    def __init__(self,comparison,unit_properties=[]):
        self._comparison=comparison
        self._unit_properties=unit_properties
        for i in range(len(self._unit_properties)):
            prop=self._unit_properties[i]
            if type(prop)!=dict:
                self._unit_properties[i]={'name':prop}
    def display(self):
        SC=self._comparison
        rows=[]
        for unit1 in SC.getSorting1().getUnitIds():
            unit2=SC.getBestUnitMatch1(unit1)
            row0={
                'Unit ID':unit1,
                'Accuracy':SC.getAgreementFraction(unit1,unit2),
                '# matches':SC.getMatchingEventCount(unit1,unit2),
                'f.n.':SC.getFalseNegativeFraction(unit1),
                'f.p.':SC.getFalsePositiveFraction(unit1),
            }
            for prop in self._unit_properties:
                pname=prop['name']
                row0[pname]=SC.getSorting1().getUnitProperty(unit_id=int(unit1),property_name=pname)
            rows.append(row0)            

        df=pd.DataFrame(rows)
        fields=['Unit ID']
        fields=fields+['Accuracy','f.n.','f.p.','# matches']
        for prop in self._unit_properties:
            pname=prop['name']
            fields.append(pname)
        df=df[fields]
        df['Accuracy'] = df['Accuracy'].map('{:,.2f}'.format)
        df['f.n.'] = df['f.n.'].map('{:,.2f}'.format)
        df['f.p.'] = df['f.p.'].map('{:,.2f}'.format)
        display(HTML(df.to_html(index=False)))