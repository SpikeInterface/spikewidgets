import pandas as pd
from IPython.display import HTML


class SortingComparisonTable():
    def __init__(self, comparison, reference=1, unit_properties=[]):
        self._comparison = comparison
        self._unit_properties = unit_properties
        self._reference = reference
        for i in range(len(self._unit_properties)):
            prop = self._unit_properties[i]
            if type(prop) != dict:
                self._unit_properties[i] = {'name': prop}

    def getDataframe(self):
        SC = self._comparison
        rows = []
        if self._reference == 1:
            for u_1, unit1 in enumerate(SC.getSorting1().getUnitIds()):
                unit2 = SC.getBestUnitMatch1(unit1)
                row0 = {
                    'Unit ID': unit1,
                    'Accuracy': SC.getAgreementFraction(unit1, unit2),
                    'Best unit': unit2,
                    'Matched unit': SC.getMappedSorting1().getMappedUnitIds(unit1),
                    '# matches': SC.getMatchingEventCount(unit1, unit2),
                    'f.n.': SC.getFalseNegativeFraction(unit1),
                    'f.p.': SC.getFalsePositiveFraction(unit1),
                }
                for prop in self._unit_properties:
                    pname = prop['name']
                    row0[pname] = SC.getSorting1().getUnitProperty(unit_id=int(unit1), property_name=pname)
                rows.append(row0)
        elif self._reference == 2:
            for u_1, unit1 in enumerate(SC.getSorting2().getUnitIds()):
                unit2 = SC.getBestUnitMatch2(unit1)
                row0 = {
                    'Unit ID': unit1,
                    'Accuracy': SC.getAgreementFraction(unit2, unit1),
                    'Best unit': unit2,
                    'Matched unit': SC.getMappedSorting2().getMappedUnitIds(unit1),
                    '# matches': SC.getMatchingEventCount(unit2, unit1),
                    'f.n.': SC.getFalseNegativeFraction(unit1),
                    'f.p.': SC.getFalsePositiveFraction(unit1),
                }
                for prop in self._unit_properties:
                    pname = prop['name']
                    row0[pname] = SC.getSorting2().getUnitProperty(unit_id=int(unit1), property_name=pname)
                rows.append(row0)

        df = pd.DataFrame(rows)
        fields = ['Unit ID']
        fields = fields + ['Accuracy', 'Best unit', 'Matched unit', 'f.n.', 'f.p.', '# matches']
        for prop in self._unit_properties:
            pname = prop['name']
            fields.append(pname)
        df = df[fields]
        df['Accuracy'] = df['Accuracy'].map('{:,.2f}'.format)
        # df['Best match'] = df['Accuracy'].map('{:,.2f}'.format)
        df['f.n.'] = df['f.n.'].map('{:,.2f}'.format)
        df['f.p.'] = df['f.p.'].map('{:,.2f}'.format)
        return df

    def display(self):
        df = self.getDataframe()
        display(HTML(df.to_html(index=False)))
