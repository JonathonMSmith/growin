import growin._sami3_utils as sami3
import growin._sami2_utils as sami2

def slice_of_growth(sami):

    import sami2py
    if isinstance(sami, sami2py._core_class.Model):
        print('Using SAMI2 model output')
        sami2.get_growth(sami)
    else:
        print('Using SAMI3 model output')
        sami3.get_growth(sami)
