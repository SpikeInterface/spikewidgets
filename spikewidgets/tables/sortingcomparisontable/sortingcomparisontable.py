import pandas as pd
from IPython.display import HTML

class SortingComparisonTable():
    def __init__(self,comparison):
        self._comparison=comparison
    def display(self):
        SC=self._comparison
        rows=[]
        for unit1 in SC.getSorting1().getUnitIds():
            unit2=SC.getBestUnitMatch1(unit1)
            rows.append({
                'Unit ID':unit1,
                'Accuracy':SC.getAgreementFraction(unit1,unit2),
                '# matches':SC.getMatchingEventCount(unit1,unit2),
                'f.n.':SC.getFalseNegativeFraction(unit1),
                'f.p.':SC.getFalsePositiveFraction(unit1),
            })
        df=pd.DataFrame(rows)
        df=df[['Unit ID','Accuracy','f.n.','f.p.','# matches']]
        df['Accuracy'] = df['Accuracy'].map('{:,.2f}'.format)
        df['f.n.'] = df['f.n.'].map('{:,.2f}'.format)
        df['f.p.'] = df['f.p.'].map('{:,.2f}'.format)
        display(HTML(df.to_html(index=False)))