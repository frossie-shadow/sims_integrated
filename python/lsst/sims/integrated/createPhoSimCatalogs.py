from __future__ import with_statement
import os
import numpy as np
from lsst.utils import getPackageDir
from lsst.afw.cameraGeom import SCIENCE
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap

from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import cached, compound
from lsst.sims.catalogs.definitions import parallelCatalogWriter

from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogPoint,
                                                          PhoSimCatalogSersic2D,
                                                          PhoSimCatalogZPoint)

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


from lsst.sims.coordUtils import getCornerRaDec

import time

__all__ = ["create_phosim_catalogs", "trim_allowed",
           "VariablePhoSimCatalogPoint", "VariablePhoSimCatalogZPoint",
           "PhoSimCatalogSersic2D_header", "StellarReferenceCatalog",
           "GalaxyReferenceCatalog"]


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
        xpix0, ypix0 = pixelCoordsFromPupilCoords(xpup, ypup, chipName='R:2,2 S:1,1', camera=_lsst_camera,
                                                  includeDistortion=False)
        return np.array([name_list, xpix, ypix, xpix0, ypix0])


class StellarReferenceCatalog(ReferenceCatalogBase, AstrometryStars, PhotometryStars, InstanceCatalog):
    pass

class GalaxyReferenceCatalog(ReferenceCatalogBase, PhotometryGalaxies, AstrometryGalaxies, InstanceCatalog):
    pass


def trim_allowed(name_list, xpix0, ypix0, magnitude, target_name, center):
    """
    Return a numpy.where() tuple indicating the objects from a list that are within
    a reasonable distance of a given chip.  "Reasonable distance" is defined such
    as to be consistent with how PhoSim's trim.cpp works: sources are accepted if
    they are within 100 pixels of the chip being considered, or if they are within

    0.1 * 2.5^(17.0 - mag)

    pixels of that bound.

    Parameters
    ----------
    name_list is a numpy array listing the names of the chips CatSim predicts
    that the sources actually land on

    xpix0, ypix0 are the TAN_PIXEL positions of the objects on the 'R:2,2 S:1,1' chip
    (i.e. the central detector)

    magnitude is a numpy array of the magNorms of the objects

    target_name is the name of the chip being considered (e.g. 'R:2,1 S:0,0')

    center is a tuple containing the TAN_PIXEL position of the center of the current
    chip relative to the 'R:2,2 S:1,1' chip.  Comparing this point to xpix0, ypix0 is
    how the method determines if an object is close and bright enough to cast scattered
    light on the detector (again, following PhoSim's prescription)

    Output
    ------
    The result of numpy.where(), indicating which of the input objects should be included
    in the input chip's InstanceCatalog
    """

    chip_radius = np.sqrt(1999.5**2 + 2035.5**2)
    distance = np.sqrt((xpix0-center[0])**2 + (ypix0-center[1])**2)
    allowed_distance = chip_radius + 1100.0 + 0.1*np.power(2.5, 17.0-magnitude)
    return np.where(np.logical_or(np.char.rfind(name_list.astype(str), target_name)>=0,
                                  distance<allowed_distance))


def _ref_cat_name_from_obs(obs, cat_dir):
    return os.path.join(cat_dir, 'phosim_%.5f_ref.txt' % obs.mjd.TAI)

def _inst_cat_name_from_obs(obs, chip_name, cat_dir):
    mangled_name = chip_name.strip().replace(':', '').replace(',', '').replace(' ', '')
    return os.path.join(cat_dir, 'phosim_%.5f_%s_cat.txt' % (obs.mjd.TAI, mangled_name))

def _write_base_pho_sim_catalogs(obs,
                                 catalog_dict={},
                                 catalog_dir=None):

    catalog_name_list = []

    ref_name = _ref_cat_name_from_obs(obs, catalog_dir)
    if os.path.exists(ref_name):
        os.unlink(ref_name)

    db_class = list(catalog_dict.keys())[0]
    db = db_class()
    ref_cat = catalog_dict[db_class][0](db, obs_metadata=obs)
    with open(ref_name, 'w') as file_handle:
        ref_cat.write_header(file_handle)

    write_header = False
    write_mode = 'a'

    for i_cat, db_class in enumerate(catalog_dict):
        db = db_class()
        print 'doing ',db.objid
        ref_cat = catalog_dict[db_class][0](db, obs_metadata=obs)
        cat_name = os.path.join(catalog_dir, 'tmp_cat_%s_%.5f_%d.txt' % (db.objid, obs.mjd.TAI, i_cat))
        if os.path.exists(cat_name):
            os.unlink(cat_name)
        catalog_name_list.append(cat_name)
        ref_cat.inst_cat_name = cat_name.split('/')[-1]
        obj_cat = catalog_dict[db_class][1](db, obs_metadata=obs)
        obj_cat.phoSimHeaderMap = {}
        local_cat_dict = {ref_name: ref_cat, cat_name: obj_cat}
        parallelCatalogWriter(local_cat_dict, chunk_size=100000,
                              write_header=write_header, write_mode=write_mode)

    return ref_name, catalog_name_list


def create_phosim_catalogs(obs_list, catalog_dir=None, db_config=None,
                           catalog_dict={StarObj: (StellarReferenceCatalog, VariablePhoSimCatalogPoint),
                                         GalaxyBulgeObj: (GalaxyReferenceCatalog, PhoSimCatalogSersic2D_header),
                                         GalaxyDiskObj: (GalaxyReferenceCatalog, PhoSimCatalogSersic2D_header),
                                         GalaxyAgnObj: (GalaxyReferenceCatalog, VariablePhoSimCatalogZPoint)}):

    t_start = time.time()

    if db_config is not None:
        config = BaseCatalogConfig()
        config.load(db_config)
        for db_class in catalog_dict:
            db_class.host = config.host
            db_class.port = config.port
            db_class.database = config.database
            db_class.driver = config.driver

    if catalog_dir is None:
        raise RuntimeError("Need to specify directory to put catalogs in")

    if not os.path.exists(catalog_dir):
        os.mkdir(catalog_dir)

    for obs in obs_list:
        ref_name, catalog_name_list = _write_base_pho_sim_catalogs(obs, catalog_dict=catalog_dict,
                                                                   catalog_dir=catalog_dir)


        print 'wrote catalog in ',time.time()-t_start

        detector_centers = {}
        for det in _lsst_camera:
            if det.getType() == SCIENCE:
                name = det.getName()
                corners = getCornerRaDec(name, _lsst_camera, obs, includeDistortion=False)

                ra_center = 0.25*(corners[0][0] + corners[1][0]
                                  + corners[2][0] + corners[3][0])

                dec_center = 0.25*(corners[0][1] + corners[1][1]
                                   + corners[2][1] + corners[3][1])

                xpix_0, ypix_0 = pixelCoordsFromRaDecLSST(ra_center, dec_center,
                                                          obs_metadata=obs,
                                                          chipName='R:2,2 S:1,1',
                                                          includeDistortion=False)

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

        star_fmt = '%s %ld %.9g %.9g %.9g %s %d %d %d %d %.9g %.9g %s %s %s %.9g %.9g\n'

        agn_dtype = np.dtype([('prefix', str, 10), ('id', long), ('ra', float), ('dec', float),
                              ('magNorm', float), ('sedFilePath', str, 200), ('redshift', float),
                              ('shear1', int), ('shear2', int), ('kappa', int), ('raOffset', float),
                              ('decOffset', float), ('spatialModel', str, 5),
                              ('internalDustModel', str, 4), ('galDustModel', str, 10), ('av', float),
                              ('rv', float)])

        agn_fmt = '%s %ld %.9g %.9g %.9g %s %.9g %d %d %d %.9g %.9g %s %s %s %.9g %.9g\n'

        gal_dtype = np.dtype([('prefix', str, 10), ('id', long), ('ra', float), ('dec', float),
                              ('magNorm', float), ('sedFilePath', str, 200), ('redshift', float),
                              ('shear1', float), ('shear2', float), ('kappa', float),
                              ('raOffset', float), ('decOffset', float), ('spatialModel', str, 8),
                              ('major', float), ('minor', float), ('pos_angle', float), ('index', float),
                              ('internalDustModel', str, 4), ('intAv', float), ('intRv', float),
                              ('galDustModel', str, 4), ('galAv', float), ('galRv', float)])

        gal_fmt = '%s %ld %.9g %.9g %.9g %s %.9g %.9g %.9g %.9g %.9g %.9g %s %.9g %.9g %.9g %.9g %s %.9g %.9g %s %.9g %.9g\n'

        db_class = list(catalog_dict.keys())[0]
        db = db_class()
        dummy_cat = catalog_dict[db_class][1](db, obs_metadata=obs)
        inst_cat_written = []

        skip_header = 1
        ct_in = 0
        catalogs_read = []
        obj_skip_dict = {}
        while skip_header == 1 or len(ref_data) == chunk_size:
            ref_data = np.genfromtxt(ref_name, dtype=ref_dtype, delimiter='; ',
                                     skip_header=skip_header, max_rows=chunk_size)
            ct_in += len(ref_data)
            skip_header += chunk_size

            temp_cat_list = np.sort(np.unique(ref_data['cat_name']))[::-1]
            for temp_cat_name in temp_cat_list:
                if 'star' in temp_cat_name.lower():
                    obj_dtype = star_dtype
                    out_fmt = star_fmt
                elif 'agn' in temp_cat_name.lower():
                    obj_dtype = agn_dtype
                    out_fmt = agn_fmt
                elif 'gal' in temp_cat_name.lower():
                    obj_dtype = gal_dtype
                    out_fmt = gal_fmt
                else:
                    raise RuntimeError('No dtype for %s' % temp_cat_name)

                valid = np.where(np.char.rfind(temp_cat_name, ref_data['cat_name'])>=0)
                local_ref_data = ref_data[valid]

                if temp_cat_name not in catalogs_read:
                    catalogs_read.append(temp_cat_name)
                    obj_skip_dict[temp_cat_name] = 0

                max_rows = len(local_ref_data)
                obj_skip = obj_skip_dict[temp_cat_name]

                obj_data = np.genfromtxt(os.path.join(catalog_dir, temp_cat_name), dtype=obj_dtype, delimiter=' ',
                                         skip_header=obj_skip, max_rows=max_rows)

                obj_skip_dict[temp_cat_name] += max_rows

                np.testing.assert_array_equal(local_ref_data['id'], obj_data['id'])

                for chip_name in detector_centers:
                    center = detector_centers[chip_name]
                    in_trim = trim_allowed(local_ref_data['chip'],
                                           local_ref_data['xpix0'], local_ref_data['ypix0'],
                                           local_ref_data['magNorm'], chip_name, center)

                    if len(in_trim[0])>0:
                        inst_cat_name = _inst_cat_name_from_obs(obs, chip_name, catalog_dir)
                        if inst_cat_name not in inst_cat_written:
                            inst_cat_written.append(inst_cat_name)
                            with open(inst_cat_name, 'w') as file_handle:
                                dummy_cat.write_header(file_handle)

                        valid_rows = obj_data[in_trim]
                        with open(inst_cat_name, 'a') as file_handle:
                            file_handle.writelines(out_fmt % tuple(vv for vv in row) for row in valid_rows)

                        print chip_name, len(in_trim[0]), len(local_ref_data), temp_cat_name.split('/')[-1]

        print ref_name,' ',ct_in

        print 'that took ',time.time()-t_start
