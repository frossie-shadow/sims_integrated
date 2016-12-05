import os
import numpy as np
import time
from lsst.utils import getPackageDir
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.photUtils import cache_LSST_seds

from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catalogs.definitions import parallelCatalogWriter

from lsst.sims.catUtils.baseCatalogModels import (BaseCatalogConfig, StarObj,
                                                  GalaxyBulgeObj, GalaxyDiskObj,
                                                  GalaxyAgnObj)

from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catUtils.mixins import AstrometryStars, AstrometryGalaxies
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies
from lsst.sims.photUtils import Sed
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST, pixelCoordsFromPupilCoords
from lsst.sims.coordUtils import _lsst_camera
from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogPoint,
                                                          PhoSimCatalogSersic2D,
                                                          PhoSimCatalogZPoint)


__all__ = ["CreatePhoSimCatalogs"]

class VariablePhoSimCatalogPoint(VariabilityStars, PhoSimCatalogPoint):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class VariablePhoSimCatalogZPoint(VariabilityStars, PhoSimCatalogZPoint):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class PhoSimCatalogSersic2D_header(PhoSimCatalogSersic2D):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class ReferenceCatalogBase(object):
    column_outputs = ['uniqueId', 'obj_type', 'raICRS', 'decICRS', 'magnitude', 'flux',
                      'chip', 'xpix', 'ypix']

    transformations = {'raICRS': np.degrees, 'decICRS':np.degrees}

    @cached
    def get_obj_type(self):
        return np.array([self.db_obj.objid]*len(self.column_by_name('raJ2000')))

    @cached
    def get_magnitude(self):
        return self.column_by_name('lsst_%s' % self.obs_metadata.bandpass)

    @cached
    def get_flux(self):
        ss = Sed()
        return ss.fluxFromMag(self.column_by_name('magnitude'))

    @compound('chip', 'xpix', 'ypix')
    def get_camera_values(self):
        xpup = self.column_by_name('x_pupil')
        ypup = self.column_by_name('y_pupil')

        name_list = chipNameFromPupilCoordsLSST(xpup, ypup)
        xpix, ypix = pixelCoordsFromPupilCoords(xpup, ypup, chipName=name_list, camera=_lsst_camera)
        return np.array([name_list, xpix, ypix])


class StellarReferenceCatalog(ReferenceCatalogBase, AstrometryStars, PhotometryStars, InstanceCatalog):
    pass

class GalaxyReferenceCatalog(ReferenceCatalogBase, PhotometryGalaxies, AstrometryGalaxies, InstanceCatalog):

    @cached
    def get_magnitude(self):
        bandpass = self.obs_metadata.bandpass
        name_lower = self.db_obj.objid.lower()
        if 'bulge' in name_lower:
            col_name = '%sBulge' % bandpass
        elif 'disk' in name_lower:
            col_name = '%sDisk' % bandpass
        elif 'agn' in name_lower:
            col_name = '%sAgn' % bandpass
        else:
            raise RuntimeError('Not sure how to get magnitudes for '
                               'db_obj: %s' % self.db_obj.objid)
        return self.column_by_name(col_name)


def CreatePhoSimCatalogs(obs_list,
                         celestial_type=('stars', 'galaxies', 'agn'),
                         catalog_dir=None):

    t_start = time.time()
    cache_LSST_seds()
    print "\ndone caching in %e\n" % (time.time()-t_start)

    config_name = os.path.join(getPackageDir('sims_integrated'), 'config', 'db.py')
    config = BaseCatalogConfig()
    config.load(config_name)
    for db_class in (StarObj, GalaxyBulgeObj, GalaxyDiskObj, GalaxyAgnObj):
        db_class.host = config.host
        db_class.port = config.port
        db_class.database = config.database
        db_class.driver = config.driver

    pkg_dir = getPackageDir('sims_integrated')
    cat_dir = os.path.join(pkg_dir, 'catalogs')
    if catalog_dir is not None:
        cat_dir = os.path.join(cat_dir, catalog_dir)
        if not os.path.exists(cat_dir):
            os.mkdir(cat_dir)

    cat_name_list = []

    for obs in obs_list:
        cat_name = os.path.join(cat_dir, 'phosim_%.5f_cat.txt' % obs.mjd.TAI)
        ref_name = os.path.join(cat_dir, 'phosim_%.5f_ref.txt' % obs.mjd.TAI)
        write_header = True
        write_mode = 'w'

        if 'stars' in celestial_type:
            db = StarObj()
            star_cat = VariablePhoSimCatalogPoint(db, obs_metadata=obs)
            ref_cat = StellarReferenceCatalog(db, obs_metadata=obs)

            cat_dict = {cat_name: star_cat, ref_name: ref_cat}

            parallelCatalogWriter(cat_dict, chunk_size=10000,
                                  write_header=write_header, write_mode=write_mode)
            write_header = False
            write_mode = 'a'

        if 'galaxies' in celestial_type:

            cat_dict = {cat_name: PhoSimCatalogSersic2D_header,
                        ref_name: GalaxyReferenceCatalog}

            for db in (GalaxyBulgeObj(), GalaxyDiskObj()):
                gal_cat = PhoSimCatalogSersic2D_header(db, obs_metadata=obs)
                ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)

                cat_dict = {cat_name: gal_cat, ref_name: ref_cat}

                parallelCatalogWriter(cat_dict, chunk_size=10000,
                                      write_header=write_header, write_mode=write_mode)

                write_header = False
                write_mode = 'a'

            db = GalaxyAgnObj()
            agn_cat = VariablePhoSimCatalogZPoint(db, obs_metadata=obs)
            ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)
            cat_dict = {cat_name: agn_cat,
                        ref_name: ref_cat}

            parallelCatalogWriter(cat_dict, chunk_size=10000,
                                  write_header=write_header, write_mode=write_mode)

        cat_name_list.append(cat_name)

    return cat_name_list
