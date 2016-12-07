from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.integrated import create_phosim_catalogs
from lsst.utils import getPackageDir

import os

import time

import copy

if __name__ == "__main__":

    t_start = time.time()

    opsimdb = os.path.join("/Users", "danielsf", "physics",
                           "lsst_150412", "Development", "garage",
                           "OpSimData", "kraken_1042_sqlite.db")

    pkg_dir = getPackageDir('sims_integrated')
    cat_dir = os.path.join(pkg_dir, 'catalogs', 'trial_161207')

    config_file = os.path.join(getPackageDir('sims_integrated'), 'config', 'db.py')

    gen = ObservationMetaDataGenerator(database=opsimdb)
    obs_list = gen.getObservationMetaData(fieldRA=(52.9, 53.1),
                                          fieldDec=(-27.5, -27.3),
                                          boundLength=0.1,
                                          telescopeFilter='r')

    print len(obs_list)

    good_obs = [obs_list[0]]

    create_phosim_catalogs(good_obs, catalog_dir=cat_dir, db_config=config_file)
    print "that took ",time.time()-t_start
