from __future__ import with_statement
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
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST, pixelCoordsFromRaDecLSST
from lsst.sims.coordUtils import pixelCoordsFromPupilCoords
from lsst.sims.coordUtils import _lsst_camera
from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogPoint,
                                                          PhoSimCatalogSersic2D,
                                                          PhoSimCatalogZPoint)


from lsst.sims.coordUtils import getCornerRaDec

import time

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
                      'chip', 'xpix', 'ypix', 'xpix0', 'ypix0', 'magNorm', 'inst_cat_name']

    transformations = {'raICRS': np.degrees, 'decICRS':np.degrees}

    delimiter = '; '

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


def _ref_cat_name_from_obs(obs, cat_dir):
    return os.path.join(cat_dir, 'phosim_%.5f_ref.txt' % obs.mjd.TAI)

def _write_base_pho_sim_catalogs(obs,
                                 celestial_type=('stars', 'galaxies', 'agn'),
                                 catalog_dir=None):

    catalog_name_dict = {}

    ref_name = _ref_cat_name_from_obs(obs, catalog_dir)
    if os.path.exists(ref_name):
        os.unlink(ref_name)

    db = StarObj()
    ref_cat = StellarReferenceCatalog(db, obs_metadata=obs)
    with open(ref_name, 'w') as file_handle:
        ref_cat.write_header(file_handle)

    catalog_name_dict[ref_name] = []
    write_header = False
    write_mode = 'a'

    if 'stars' in celestial_type:
        db = StarObj()
        ref_cat = StellarReferenceCatalog(db, obs_metadata=obs)
        cat_name = os.path.join(catalog_dir, 'tmp_stars_phosim_%.5f_cat.txt' % obs.mjd.TAI)
        if os.path.exists(cat_name):
            os.unlink(cat_name)
        catalog_name_dict[ref_name].append(cat_name)
        ref_cat.inst_cat_name = cat_name.split('/')[-1]
        star_cat = VariablePhoSimCatalogPoint(db, obs_metadata=obs)
        star_cat.phoSimHeaderMap = {}
        cat_dict = {ref_name: ref_cat, cat_name: star_cat}

        parallelCatalogWriter(cat_dict, chunk_size=100000,
                              write_header=write_header, write_mode=write_mode)
        write_mode = 'a'
        write_header = False
        print 'done with stars'

    if 'galaxies' in celestial_type:

        cat_name = os.path.join(catalog_dir, 'tmp_galaxies_phosim_%.5f_cat.txt' % obs.mjd.TAI)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        catalog_name_dict[ref_name].append(cat_name)

        for db in (GalaxyBulgeObj(), GalaxyDiskObj()):
            ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)

            ref_cat.inst_cat_name = cat_name.split('/')[-1]
            gal_cat = PhoSimCatalogSersic2D_header(db, obs_metadata=obs)
            gal_cat.phoSimHeaderMap = {}
            cat_dict = {ref_name: ref_cat, cat_name: gal_cat}

            parallelCatalogWriter(cat_dict, chunk_size=100000,
                                  write_header=write_header, write_mode=write_mode)

            write_header = False
            write_mode = 'a'
            print 'done with ',db.objid

        db = GalaxyAgnObj()
        ref_cat = GalaxyReferenceCatalog(db, obs_metadata=obs)
        cat_name = os.path.join(catalog_dir, 'tmp_agn_phosim_%.5f_cat.txt' % obs.mjd.TAI)
        if os.path.exists(cat_name):
            os.unlink(cat_name)
        catalog_name_dict[ref_name].append(cat_name)
        ref_cat.inst_cat_name = cat_name.split('/')[-1]
        gal_cat = VariablePhoSimCatalogZPoint(db, obs_metadata=obs)
        gal_cat.phoSimHeaderMap = {}
        cat_dict = {ref_name: ref_cat, cat_name: gal_cat}

        parallelCatalogWriter(cat_dict, chunk_size=100000,
                              write_header=write_header, write_mode=write_mode)

    return catalog_name_dict


def CreatePhoSimCatalogs(obs_list, celestial_type=('stars', 'galaxies', 'agn'),
                         catalog_dir=None):

    t_start = time.time()

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

    for obs in obs_list:
        catalog_name_dict = _write_base_pho_sim_catalogs(obs, celestial_type=celestial_type,
                                                         catalog_dir=cat_dir)


        detector_centers = {}
        for det in _lsst_camera:
            if det.getType() == SCIENCE:
                name = det.getName()
                corners = getCornerRaDec(name, _lsst_camera, obs)

                ra_center = 0.25*(corners[0][0] + corners[1][0]
                                  + corners[2][0] + corners[3][0])

                dec_center = 0.25*(corners[0][1] + corners[1][1]
                                   + corners[2][1] + corners[3][1])

                xpix_0, ypix_0 = pixelCoordsFromRaDecLSST(ra_center, dec_center,
                                                          obs_metadata=obs,
                                                          chipName='R:2,2 S:1,1')

                detector_centers[name] = (xpix_0, ypix_0)

        # code to read in the reference catalog
        chunk_size = 100000
        ref_dtype = np.dtype([('id', long), ('type', str, 100), ('ra', float), ('dec', float),
                              ('chip', str, 100), ('xpix', float), ('ypix', float),
                              ('xpix0', float), ('ypix0', float), ('magNorm', float),
                              ('cat_name', str, 100)])

        star_dtype = np.dtype([('prefix', str, 10), ('id', long), ('ra', float), ('dec', float),
                               ('magNorm', float), ('sedFilePath', str, 200), ('redshift', int),
                               ('shear1', int), ('shear2', int), ('kappa', int), ('raOffset', float),
                               ('decOffset', float), ('spatialModel', str, 5),
                               ('internalDustModel', str, 4), ('galDustModel', str, 10), ('av', float),
                               ('rv', float)])

        agn_dtype = np.dtype([('prefix', str, 10), ('id', long), ('ra', float), ('dec', float),
                              ('magNorm', float), ('sedFilePath', str, 200), ('redshift', float),
                              ('shear1', int), ('shear2', int), ('kappa', int), ('raOffset', float),
                              ('decOffset', float), ('spatialModel', str, 5),
                              ('internalDustModel', str, 4), ('galDustModel', str, 10), ('av', float),
                              ('rv', float)])

        gal_dtype = np.dtype([('prefix', str, 10), ('id', long), ('ra', float), ('dec', float),
                              ('magNorm', float), ('sedFilePath', str, 200), ('redshift', float),
                              ('shear1', float), ('shear2', float), ('kappa', float),
                              ('raOffset', float), ('decOffset', float), ('spatialModel', str, 8),
                              ('internalDustModel', str, 4), ('intAv', float), ('intRv', float),
                              ('galDustModel', str, 4), ('galAv', float), ('galRv', float)])


        for ref_cat in catalog_name_dict:
            skip_header = 1
            ct_in = 0
            catalogs_written = []
            obj_skip_dict = {}
            while skip_header == 1 or len(ref_data) == chunk_size:
                ref_data = np.genfromtxt(ref_cat, dtype=ref_dtype, delimiter='; ',
                                         skip_header=skip_header, max_rows=chunk_size)
                ct_in += len(ref_data)
                skip_header += chunk_size

                inst_cat_list = np.unique(ref_data['cat_name'])
                for cat_name in inst_cat_list:
                    if 'star' in cat_name:
                        obj_dtype = star_dtype
                    elif 'agn' in cat_name:
                        obj_dtype = agn_dtype
                    elif 'gal' in cat_name:
                        obj_dtype = gal_dtype
                    else:
                        raise RuntimeError('No dtype for %s' % cat_name)

                    valid = np.where(np.char.rfind(cat_name, ref_data['cat_name'])>=0)
                    local_ref_data = ref_data[valid]

                    if cat_name not in catalogs_written:
                        catalogs_written.append(cat_name)
                        obj_skip_dict[cat_name] = 0

                    max_rows = len(local_ref_data)
                    obj_skip = obj_skip_dict[cat_name]

                    obj_data = np.genfromtxt(os.path.join(cat_dir, cat_name), dtype=obj_dtype, delimiter=' ',
                                             skip_header=obj_skip, max_rows=max_rows)

                    obj_skip_dict[cat_name] += max_rows

                    np.testing.assert_array_equal(local_ref_data['id'], obj_data['id'])


            print ref_cat,' ',ct_in

        print 'that took ',time.time()-t_start
