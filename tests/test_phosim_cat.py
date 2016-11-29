from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.integrated import CreatePhoSimCatalogs

import os

import time

import copy

if __name__ == "__main__":

    t_start = time.time()

    opsimdb = os.path.join("/Users", "danielsf", "physics",
                           "lsst_150412", "Development", "garage",
                           "OpSimData", "kraken_1042_sqlite.db")

    gen = ObservationMetaDataGenerator(database=opsimdb)
    obs_list = gen.getObservationMetaData(fieldRA=(52.9, 53.1),
                                          fieldDec=(-27.5, -27.3),
                                          boundLength=0.5,
                                          telescopeFilter='r')

    print len(obs_list)

    good_obs = [obs_list[0], obs_list[-1]]

    CreatePhoSimCatalogs(good_obs, catalog_dir='trial')
    print "that took ",time.time()-t_start
