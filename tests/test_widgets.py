import spikeextractors as se
import spiketoolkit as st
import spikewidgets as sw
import spikecomparison
import unittest


class TestWidgets(unittest.TestCase):
    def setUp(self):
        self._RX, self._SX = se.example_datasets.toy_example(num_channels=4, duration=10)

    def tearDown(self):
        pass

    def test_timeseries(self):
        sw.plot_timeseries(self._RX)

    def test_geometry(self):
        sw.plot_electrode_geometry(self._RX)

    def test_unitwaveforms(self):
        sw.plot_unit_waveforms(self._RX, self._SX)

    def test_ampdist(self):
        sw.plot_amplitudes_distribution(self._RX, self._SX)

    def test_amptime(self):
        sw.plot_amplitudes_timeseries(self._RX, self._SX)

    def test_features(self):
        sw.plot_features(self._RX, self._SX)

    def test_ach(self):
        sw.plot_autocorrelograms(self._SX, bin_size=1, window=10)

    def test_cch(self):
        sw.plot_crosscorrelograms(self._SX, bin_size=1, window=10)

    def test_isi(self):
        sw.plot_isi_distribution(self._SX, bins=10, window=1)

    def test_rasters(self):
        sw.plot_rasters(self._SX)

    def test_confusion(self):
        sc = spikecomparison.compare_two_sorters(self._SX, self._SX)
        sw.plot_confusion_matrix(sc, count_text=False)

    def test_multicomp_graph(self):
        msc = spikecomparison.compare_multiple_sorters([self._SX, self._SX, self._SX])
        sw.plot_multicomp_graph(msc, edge_cmap='viridis', node_cmap='rainbow', draw_labels=False)


if __name__ == '__main__':
    unittest.main()
