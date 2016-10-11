from lsst.sims.catUtils.baseCatalogModels import (StarObj,
                                                  GalaxyBulgeObj, GalaxyDiskObj,
                                                  GalaxyAgnObj)

from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogPoint,
                                                          PhoSimCatalogSersic2D,
                                                          PhoSimCatalogZPoint)

__all__ = ["celestial_db_dict"]

class VariablePhoSimCatalogPoint(VariabilityStars, PhoSimCatalogPoint):
    pass

class VariablePhoSimCatalogZPoint(VariabilityStars, PhoSimCatalogZPoint):
    pass

celestial_db_dict = {'stars': ([StarObj], [VariablePhoSimCatalogPoint]),
                     'galaxies': ([GalaxyBulgeObj, GalaxyDiskObj],
                                  [PhoSimCatalogSersic2D, PhoSimCatalogSersic2D]),
                     'agn': ([GalaxyAgnObj], [VariablePhoSimCatalogZPoint])}
