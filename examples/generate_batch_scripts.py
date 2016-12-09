import os
from lsst.utils import getPackageDir
from lsst.sims.integrated import create_phosim_catalogs, create_bash_scripts
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

opsimdb = os.path.join("/Users", "danielsf", "physics",
                       "lsst_150412", "Development", "garage",
                       "OpSimData", "kraken_1042_sqlite.db")

gen = ObservationMetaDataGenerator(database=opsimdb)
obs_list = gen.getObservationMetaData(fieldRA=(52.9, 53.1),
                                      fieldDec=(-27.5, -27.3),
                                      boundLength=2.0,
                                      telescopeFilter='r')


db_config = os.path.join(getPackageDir('sims_integrated'), 'config',
                         'db.py')

catalog_dir = os.path.join(getPackageDir('sims_integrated'), 'examples', 'cat')

create_phosim_catalogs([obs_list[0]], catalog_dir=catalog_dir, db_config=db_config)

output_dir = os.path.join(getPackageDir('sims_integrated'), 'examples', 'images')

batch_dir = os.path.join(getPackageDir('sims_integrated'), 'examples', 'batch_dir')

create_bash_scripts(catalog_dir, n_per_batch=189/40,
                    output_dir=output_dir, batch_dir=batch_dir,
                    db_config=db_config)
