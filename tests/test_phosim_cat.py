from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.integrated import CreatePhoSimCatalogs

import os

import time

if __name__ == "__main__":

    t_start = time.time()

    opsimdb = os.path.join("/Users", "danielsf", "physics",
                           "lsst_150412", "Development", "garage",
                           "OpSimData", "kraken_1042_sqlite.db")

    gen = ObservationMetaDataGenerator(database=opsimdb)
    obs_list = gen.getObservationMetaData(fieldRA=(23.0, 50.0),
                                          fieldDec=(-43.0, -39.0),
                                          expMJD=(59580.0, 60000.0))

    print len(obs_list)

    CreatePhoSimCatalogs(obs_list[:2], catalog_dir='trial')
    print "that took ",time.time()-t_start
