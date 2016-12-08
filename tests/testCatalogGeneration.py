from __future__ import with_statement
import unittest
import os
import numpy as np

import lsst.utils.tests
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.utils import getPackageDir
from lsst.sims.integrated import trim_allowed, create_phosim_catalogs
from lsst.sims.catalogs.definitions import CompoundInstanceCatalog
from lsst.sims.utils import ObservationMetaData
from lsst.sims.coordUtils import pixelCoordsFromPupilCoords, getCornerRaDec, _lsst_camera
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.coordUtils import pixelCoordsFromRaDecLSST
from lsst.sims.catalogs.decorators import cached
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.utils import makePhoSimTestDB, TestVariabilityMixin
from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogPoint
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogZPoint
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catUtils.baseCatalogModels import (StarObj, GalaxyBulgeObj, GalaxyDiskObj,
                                                  GalaxyAgnObj)

from lsst.sims.integrated import (StellarReferenceCatalog, GalaxyReferenceCatalog,
                                  VariablePhoSimCatalogPoint, VariablePhoSimCatalogZPoint,
                                  PhoSimCatalogSersic2D_header)

class TestDBObj(object):

    def query_columns(*args, **kwargs):
        return CatalogDBObject.query_columns(*args, **kwargs)

    def _get_column_query(*args, **kwargs):
        return CatalogDBObject._get_column_query(*args, **kwargs)
    
    def _final_pass(*args, **kwargs):
        return CatalogDBObject._final_pass(*args, **kwargs)


class TestStarObj(TestDBObj, StarObj, CatalogDBObject):
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

class TestGalaxyBulgeObj(TestDBObj, GalaxyBulgeObj, CatalogDBObject):
    pass

class TestGalaxyDiskObj(TestDBObj, GalaxyDiskObj, CatalogDBObject):
    pass

class TestAgnObj(TestDBObj, GalaxyAgnObj, CatalogDBObject):
    pass

class TestPhoSimPoint(TestVariabilityMixin, VariablePhoSimCatalogPoint):
    pass

class TestPhoSimZPoint(TestVariabilityMixin, VariablePhoSimCatalogZPoint):
    pass

class TestPhoSimSersic2D(TestVariabilityMixin, PhoSimCatalogSersic2D_header):
    pass   

class PhoSimControl(object):
    phoSimHeaderMap = DefaultPhoSimHeaderMap
    chip_name = None
    cannot_be_null = ['sedFilepath', 'onChip']

    @cached
    def get_onChip(self):

        xpup = self.column_by_name('x_pupil')
        ypup = self.column_by_name('y_pupil')
        magNorm = self.column_by_name('phoSimMagNorm')
        if len(xpup) == 0:
            return np.array([])

        if not hasattr(self, 'chip_center'):
            corners = getCornerRaDec(self.chip_name, _lsst_camera, self.obs_metadata,
                                     includeDistortion=False)
            ra_center = 0.25*(corners[0][0] + corners[1][0]
                              + corners[2][0] + corners[3][0])

            dec_center = 0.25*(corners[0][1] + corners[1][1]
                               + corners[2][1] + corners[3][1])

            self.chip_center = pixelCoordsFromRaDecLSST(ra_center, dec_center,
                                                        obs_metadata=self.obs_metadata,
                                                        chipName='R:2,2 S:1,1',
                                                        includeDistortion=False)

        xpix, ypix = pixelCoordsFromPupilCoords(xpup, ypup, chipName='R:2,2 S:1,1',
                                                camera=_lsst_camera, includeDistortion=False)

        name_list = chipNameFromPupilCoordsLSST(xpup, ypup)

        dexes = trim_allowed(name_list, xpix, ypix, magNorm, self.chip_name, self.chip_center)
        return np.array(['allowed' if ii in dexes[0] else 'NULL' for ii in range(len(xpup))])


class PhoSimStarControl(PhoSimControl, TestVariabilityMixin, VariabilityStars, PhoSimCatalogPoint):
    pass


class PhoSimAgnControl(PhoSimControl, TestVariabilityMixin, VariabilityStars, PhoSimCatalogZPoint):
    pass


class PhoSimGalControl(PhoSimControl, TestVariabilityMixin, PhoSimCatalogSersic2D):
    pass


def setup_module(module):
    lsst.utils.tests.init()


class PhoSimCatalogCreationTestCase(unittest.TestCase):

    longMessage = True

    @classmethod
    def setUpClass(cls):
        pkg_dir = getPackageDir('sims_integrated')
        cls.scratch_dir = os.path.join(pkg_dir, 'tests', 'scratchSpace')

        cls.db_name = os.path.join(cls.scratch_dir, 'sims_integrated_test_db.db')
        if os.path.exists(cls.db_name):
            os.unlink(cls.db_name)

        cls.obs = makePhoSimTestDB(filename=cls.db_name, size=10000, seedVal=8123361,
                                   radius = 2.0)

        cls.config_name = os.path.join(cls.scratch_dir, 'sims_integrated_test_config.py')
        with open(cls.config_name, 'w') as file_handle:
            file_handle.write("config.driver='sqlite'\n")
            file_handle.write("config.host='None'\n")
            file_handle.write("config.port='None'\n")
            file_handle.write("config.database='%s'\n" % cls.db_name)

    #@classmethod
    #def tearDownClass(cls):
    #    sims_clean_up()
    #    if os.path.exists(cls.db_name):
    #        os.unlink(cls.db_name)
    #    if os.path.exists(cls.config_name):
    #        os.unlink(cls.config_name)

    def test_catalog_generation(self):
        """
        Generate PhoSimInstance catalogs using the create_phosim_catalogs method.
        Verify that the results are the same as if we had generated the same catalogs
        looping over each chip by hand.
        """

        catalog_dir = os.path.join(self.scratch_dir, 'test_cat_gen_dir')

        catalog_dict = {TestStarObj: (StellarReferenceCatalog, TestPhoSimPoint),
                        TestGalaxyBulgeObj: (GalaxyReferenceCatalog, TestPhoSimSersic2D),
                        TestGalaxyDiskObj: (GalaxyReferenceCatalog, TestPhoSimSersic2D),
                        TestAgnObj: (GalaxyReferenceCatalog, TestPhoSimZPoint)}

        create_phosim_catalogs([self.obs], catalog_dir=catalog_dir, db_config=self.config_name,
                               catalog_dict=catalog_dict)

        for det in _lsst_camera:
            chip_name = det.getName()
            mangled_name = chip_name.strip().replace(':', '').replace(',', '').replace(' ', '')
            test_cat_name = os.path.join(catalog_dir,
                                         'phosim_%.5f_%s_cat.txt' % (self.obs.mjd.TAI,
                                                                     mangled_name))

            control_cat = CompoundInstanceCatalog([PhoSimStarControl,
                                                   PhoSimGalControl,
                                                   PhoSimGalControl,
                                                   PhoSimAgnControl],
                                                  [TestStarObj, TestGalaxyBulgeObj,
                                                   TestGalaxyDiskObj, TestAgnObj],
                                                  obs_metadata=self.obs)

            control_cat.phoSimHeaderMap = DefaultPhoSimHeaderMap
            control_cat.chip_name = chip_name

            control_cat_name = os.path.join(self.scratch_dir,
                                            "phosim_creation_control_catalog.txt")

            if os.path.exists(control_cat_name):
                os.unlink(control_cat_name)

            control_cat.write_catalog(control_cat_name)

            with open(test_cat_name, 'r') as file_handle:
                test_lines = file_handle.readlines()

            with open(control_cat_name, 'r') as file_handle:
                control_lines = file_handle.readlines()

            for line in test_lines:
                self.assertIn(line, control_lines, msg='%s\n not in control' % line)

            for line in control_lines:
                self.assertIn(line, test_lines, msg='%s\n not in test' % line)

            if os.path.exists(control_cat_name):
                os.unlink(control_cat_name)

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
