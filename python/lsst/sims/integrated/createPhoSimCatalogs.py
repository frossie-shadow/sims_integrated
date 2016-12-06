import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.afw.cameraGeom import SCIENCE
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap

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

class PhoSimTrimBase(object):

    cannot_be_null = ['sedFilepath', 'trim_allowed']
    chip_name = None

    @compound('chip', 'xpix', 'ypix')
    def get_camera_values(self):
        xpup = self.column_by_name('x_pupil')
        ypup = self.column_by_name('y_pupil')

        name_list = chipNameFromPupilCoordsLSST(xpup, ypup)
        xpix, ypix = pixelCoordsFromPupilCoords(xpup, ypup, chipName=name_list, camera=_lsst_camera)
        return np.array([name_list, xpix, ypix])

    @cached
    def get_trim_allowed(self):
        """
        Return 'allowed' for any objects predicted to be either on the current chip
        or within 100 + 0.1*2.5^(17-magNorm) pixels of the current chip (this is
        the buffer applied by PhoSim's trim.cpp)
        """
        name_list = self.column_by_name('chip')
        xpup = self.column_by_name('x_pupil')
        ypup = self.column_by_name('y_pupil')
        mag_list = self.column_by_name('magNorm')

        if len(name_list) == 0:
            return np.array([])

        if self.chip_name is None:
            raise RuntimeError("Cannot perform trimming of InstanceCatalogs; "
                               "you have not set chip_name in one of your catalogs: %s " % self.db_obj.objid)

        xpix, ypix = pixelCoordsFromPupilCoords(xpup, ypup, chipName=self.chip_name, camera=_lsst_camera)

        chip_radius = np.sqrt(1999.5**2 + 2035.5**2)

        distance = np.sqrt((xpix-1999.5)**2 + (ypix-2035.5)**2)
        allowed_distance = chip_radius + 100.0 + 0.1*np.power(2.5, 17.0-mag_list)
        return np.where(np.logical_or(np.char.rfind(name_list.astype(str), self.chip_name)>=0,
                                      distance<allowed_distance), 'valid', 'NULL')


class VariablePhoSimCatalogPoint(VariabilityStars, PhoSimCatalogPoint):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class VariablePhoSimCatalogZPoint(VariabilityStars, PhoSimCatalogZPoint):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class PhoSimCatalogSersic2D_header(PhoSimCatalogSersic2D):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class ReferenceCatalogBase(object):
    column_outputs = ['uniqueId', 'obj_type', 'raICRS', 'decICRS',
                      'chip', 'xpix', 'ypix', 'xpix0', 'ypix0', 'inst_cat_name']

    transformations = {'raICRS': np.degrees, 'decICRS':np.degrees}

    inst_cat_name = None

    @cached
    def get_obj_type(self):
        return np.array([self.db_obj.objid]*len(self.column_by_name('raJ2000')))

    @cached
    def get_inst_cat_name(self):
        return np.array([self.inst_cat_name]*len(self.column_by_name('raJ2000')))

    @compound('chip', 'xpix', 'ypix', 'xpix0', 'ypix0')
    def get_camera_values(self):
        xpup = self.column_by_name('x_pupil')
        ypup = self.column_by_name('y_pupil')

        name_list = chipNameFromPupilCoordsLSST(xpup, ypup)
        xpix, ypix = pixelCoordsFromPupilCoords(xpup, ypup, chipName=name_list, camera=_lsst_camera)
        xpix0, ypix0 = pixelCoordsFromPupilCoords(xpup, ypup, chipName='R:2,2 S:1,1', camera=_lsst_camera)
        return np.array([name_list, xpix, ypix, xpix0, ypix0])


class StellarReferenceCatalog(ReferenceCatalogBase, AstrometryStars, PhotometryStars, InstanceCatalog):
    pass

class GalaxyReferenceCatalog(ReferenceCatalogBase, PhotometryGalaxies, AstrometryGalaxies, InstanceCatalog):
    pass


def CreatePhoSimCatalogs(obs_list,
                         celestial_type=('stars', 'galaxies', 'agn'),
                         catalog_dir=None):

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
        ref_name = os.path.join(cat_dir, 'phosim_%.5f_ref.txt' % obs.mjd.TAI)
        write_header = True
        write_mode = 'w'

        if 'stars' in celestial_type:
            db = StarObj()
            ref_cat = StellarReferenceCatalog(db, obs_metadata=obs)
            cat_name = os.path.join(cat_dir, 'tmp_stars_phosim_%.5f_cat.txt' % obs.mjd.TAI)
            ref_cat.inst_cat_name = cat_name.split('/')[-1]
            star_cat = VariablePhoSimCatalogPoint(db, obs_metadata=obs)
            cat_dict = {ref_name: ref_cat, cat_name: star_cat}

            parallelCatalogWriter(cat_dict, chunk_size=10000,
                                  write_header=write_header, write_mode=write_mode)
            write_header = False
            write_mode = 'a'
            print 'done with stars'

        if 'galaxies' in celestial_type:

            for db in (GalaxyBulgeObj(), GalaxyDiskObj()):
                ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)
                cat_name = os.path.join(cat_dir, 'tmp_galaxies_phosim_%.5f_cat.txt' % obs.mjd.TAI)
                ref_cat.inst_cat_name = cat_name.split('/')[-1]
                gal_cat = PhoSimCatalogSersic2D_header(db, obs_metadata=obs)
                cat_dict = {ref_name: ref_cat, cat_name: gal_cat}

                parallelCatalogWriter(cat_dict, chunk_size=10000,
                                      write_header=write_header, write_mode=write_mode)

                write_header = False
                write_mode = 'a'
                print 'done with ',db.objid

            db = GalaxyAgnObj()
            ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)
            cat_name = os.path.join(cat_dir, 'tmp_agn_phosim_%.5f_cat.txt' % obs.mjd.TAI)
            ref_cat.inst_cat_name = cat_name.split('/')[-1]
            gal_cat = VariablePhoSimCatalogZPoint(db, obs_metadata=obs)
            cat_dict = {ref_name: ref_cat, cat_name: gal_cat}

            parallelCatalogWriter(cat_dict, chunk_size=10000,
                                  write_header=write_header, write_mode=write_mode)

        cat_name_list.append(cat_name)

    return cat_name_list
