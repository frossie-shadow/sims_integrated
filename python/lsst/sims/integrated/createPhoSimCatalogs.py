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


class VariablePhoSimCatalogPoint(VariabilityStars, PhoSimTrimBase, PhoSimCatalogPoint):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class VariablePhoSimCatalogZPoint(VariabilityStars, PhoSimTrimBase, PhoSimCatalogZPoint):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class PhoSimCatalogSersic2D_header(PhoSimTrimBase, PhoSimCatalogSersic2D):
    phoSimHeaderMap = DefaultPhoSimHeaderMap

class ReferenceCatalogBase(object):
    column_outputs = ['uniqueId', 'obj_type', 'raICRS', 'decICRS',
                      'chip', 'xpix', 'ypix', 'xpix0', 'ypix0']

    transformations = {'raICRS': np.degrees, 'decICRS':np.degrees}

    @cached
    def get_obj_type(self):
        return np.array([self.db_obj.objid]*len(self.column_by_name('raJ2000')))

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


def _append_cat_dict(cat_class, db, obs, cat_dir, cat_dict):

    if not hasattr(_append_cat_dict, 'camera_map'):
        _append_cat_dict.camera_map = {}
        for det in _lsst_camera:
            if det.getType() == SCIENCE:
                name = det.getName()
                cat_name = name.strip().replace(':', '_').replace(',', '_').replace(' ','_')
                _append_cat_dict.camera_map[name] = cat_name

    for chip_name in _append_cat_dict.camera_map:
        cat = cat_class(db, obs_metadata=obs)
        cat.chip_name = chip_name
        safe_name = _append_cat_dict.camera_map[chip_name]
        file_name = os.path.join(cat_dir,
                                 'phosim_%s_%.5f.txt' % (safe_name, obs.mjd.TAI))

        cat_dict[file_name] = cat

    return cat_dict

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
        cat_name = os.path.join(cat_dir, 'phosim_%.5f_cat.txt' % obs.mjd.TAI)
        ref_name = os.path.join(cat_dir, 'phosim_%.5f_ref.txt' % obs.mjd.TAI)
        write_header = True
        write_mode = 'w'

        if 'stars' in celestial_type:
            db = StarObj()
            ref_cat = StellarReferenceCatalog(db, obs_metadata=obs)
            cat_dict = {ref_name: ref_cat}
            cat_dict = _append_cat_dict(VariablePhoSimCatalogPoint, db, obs,
                                        cat_dir, cat_dict)

            parallelCatalogWriter(cat_dict, chunk_size=10000,
                                  write_header=write_header, write_mode=write_mode)
            write_header = False
            write_mode = 'a'
            print 'done with stars'

        if 'galaxies' in celestial_type:

            for db in (GalaxyBulgeObj(), GalaxyDiskObj()):
                ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)
                cat_dict = {ref_name: ref_cat}
                cat_dict = _append_cat_dict(PhoSimCatalogSersic2D_header, db, obs,
                                            cat_dir, cat_dict)

                parallelCatalogWriter(cat_dict, chunk_size=10000,
                                      write_header=write_header, write_mode=write_mode)

                write_header = False
                write_mode = 'a'
                print 'done with ',db.objid

            db = GalaxyAgnObj()
            ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)
            cat_dict = {ref_name: ref_cat}
            cat_dict = _append_cat_dict(VariablePhoSimCatalogZPoint, db, obs,
                                        cat_dir, cat_dict)

            parallelCatalogWriter(cat_dict, chunk_size=10000,
                                  write_header=write_header, write_mode=write_mode)

        cat_name_list.append(cat_name)

    return cat_name_list
