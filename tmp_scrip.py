import pickle
import irfl


sami = pickle.load(open('/home/jmsmit37/Projects/sami2_rtgr/data_out/growth/lon022/2011_355/sami.p', 'rb'))
irfl.generate_plots.plot_growth_term(sami, 'DecSol', drift=True)
