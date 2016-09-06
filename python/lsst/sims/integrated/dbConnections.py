from lsst.sims.catUtils.baseCatalogModels import (StarObj,
                                                  GalaxyBulgeObj, GalaxyDiskObj,
                                                  GalaxyAgnObj)

from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogPoint,
                                                          PhoSimCatalogSersic2D,
                                                          PhoSimCatalogZPoint)

__all__ = ["celestial_db_dict"]

celestial_db_dict = {'stars': ([StarObj], [PhoSimCatalogPoint]),
                     'galaxies': ([GalaxyBulgeObj, GalaxyDiskObj],
                                  [PhoSimCatalogSersic2D, PhoSimCatalogSersic2D]),
                     'agn': ([GalaxyAgnObj], [PhoSimCatalogZPoint])}
