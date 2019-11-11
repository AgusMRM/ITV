import numpy as np
import yt
from yt.analysis_modules.halo_analysis.api import HaloCatalog
import pynbody

data_ds = yt.load('/mnt/is2/dpaz/ITV/S1373/out/snapshot_050')
hc = HaloCatalog(data_ds=data_ds, finder_method='fof',
		finder_kwargs={"ptype": "stars","padding": 0.02})
hc.create()
