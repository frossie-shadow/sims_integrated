import os
from lsst.utils import getPackageDir
from lsst.sims.catUtils.baseCatalogModels import GalaxyTileCompoundObj
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.catalogs.definitions import CompoundInstanceCatalog

from lsst.sims.catUtils.baseCatalogModels import (StarObj,
                                                  GalaxyBulgeObj, GalaxyDiskObj,
                                                  GalaxyAgnObj)

from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogPoint,
                                                          PhoSimCatalogSersic2D,
                                                          PhoSimCatalogZPoint)


__all__ = ["CreatePhoSimCatalogs"]

class VariablePhoSimCatalogPoint(VariabilityStars, PhoSimCatalogPoint):
    pass

class VariablePhoSimCatalogZPoint(VariabilityStars, PhoSimCatalogZPoint):
    pass

celestial_db_dict = {'stars': ([StarObj], [VariablePhoSimCatalogPoint]),
                     'galaxies': ([GalaxyBulgeObj, GalaxyDiskObj],
                                  [PhoSimCatalogSersic2D, PhoSimCatalogSersic2D]),
                     'agn': ([GalaxyAgnObj], [VariablePhoSimCatalogZPoint])}

def CreatePhoSimCatalogs(obs_list,
                         celestial_type=('stars', 'galaxies', 'agn'),
                         catalog_dir=None):

    db_class_list = []
    cat_class_list = []
    for cc in celestial_type:
        db_class_list += celestial_db_dict[cc][0]
        cat_class_list += celestial_db_dict[cc][1]

    pkg_dir = getPackageDir('sims_integrated')
    cat_dir = os.path.join(pkg_dir, 'catalogs')
    if catalog_dir is not None:
        cat_dir = os.path.join(cat_dir, catalog_dir)
        if not os.path.exists(cat_dir):
            os.mkdir(cat_dir)

    cat_name_list = []
    compound_cat = None
    connection_list = None
    for obs in obs_list:
        if compound_cat is not None:
            connection_list = compound_cat._active_connections

        cat_name = 'phosim_%.5f_cat.txt' % obs.mjd.TAI

        compound_cat = CompoundInstanceCatalog(cat_class_list,
                                               db_class_list,
                                               obs_metadata=obs,
                                               compoundDBclass=GalaxyTileCompoundObj)

        if connection_list is not None:
            compound_cat._active_connections += connection_list

        compound_cat.phoSimHeaderMap = DefaultPhoSimHeaderMap

        compound_cat.write_catalog(os.path.join(cat_dir, cat_name), chunk_size=1000000)
        cat_name_list.append(cat_name)

    return cat_name_list
