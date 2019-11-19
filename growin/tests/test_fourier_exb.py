"""Unit testing for the functions and objects in fourier_exb
"""
import numpy as np
import pytest
import growin

class TestDriftInstrument():
    """tests of the drift instrument class
    """
    def test_model_instantiation(self):
        """instantiate a Drift
        """
        drift_inst = growin.DriftInstrument(platform='cnofs', name='ivm')
        assert isinstance(drift_inst, growin.fourier_exb.DriftInstrument)

    def test_model_pipeline(self):
        """calculate the median drifts for a single region and check output
        """
        drift_inst = growin.DriftInstrument(platform='cnofs', name='ivm',
                                            clean_level='none',
                                            strict_time_flag=True)
        drift_inst.get_drifts(drift_key='ionVelmeridional', start_year=2010,
                              stop_year=2010)
        assert hasattr(drift_inst, 'drifts')
        assert hasattr(drift_inst.drifts, 'drifts')
        assert isinstance(drift_inst.drifts.drifts.values, (np.ndarray))
        assert drift_inst.drifts.drifts.values.shape == (1, 1, 1, 48)


    def test_drift_key_error(self):
        """if there is no drift key supplied to get_drifts it should raise
           an error
        """
        with pytest.raises(AttributeError):
            drift_inst = growin.DriftInstrument(platform='cnofs', name='ivm')
            drift_inst.get_drifts()
