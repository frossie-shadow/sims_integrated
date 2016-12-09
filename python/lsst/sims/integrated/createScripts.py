import os
import re
from lsst.utils import getPackageDir


def create_bash_scripts(catalog_dir, n_per_batch, output_dir, batch_dir,
                        physics_file=None, db_config=None):

    if not os.path.exists(batch_dir):
        os.mkdir(batch_dir)

    sed_dir = getPackageDir('sims_sed_library')

    inst_cat_list = os.listdir(catalog_dir)

    bash_name_list = []

    n_cmd = 0
    n_file = 0
    file_handle = None
    for inst_cat in inst_cat_list:
        if inst_cat.endswith('_cat.txt'):
            if file_handle is None:
                file_handle = open(os.path.join(batch_dir, "phosim_bash_%d.sh" % n_file), 'w')
                n_file += 1

            raw_name = re.search('R[0-9][0-9]S[0-9][0-9]', inst_cat).group(0)
            chip_name = '%s:%s,%s %s:%s,%s' % tuple(rr for rr in raw_name)

            file_handle.write('./phosim %s' % os.path.join(catalog_dir, inst_cat))
            file_handle.write(' -o %s' % output_dir)
            file_handle.write(' -sed=%s' % sed_dir)
            file_handle.write(' -s %s' % chip_name)
            if physics_file is not None:
                file_handle.write(' -c %s' % physics_file)
            file_handle.write('\n\n')

            n_cmd += 1
            if n_cmd == n_per_batch:
                file_handle.close()
                file_handle = None
                n_cmd = 0

    if file_handle is not None:
        file_handle.close()
