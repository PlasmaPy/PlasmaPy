import os
import glob

import plasmapy

rootdir = os.path.join(os.path.dirname(plasmapy.__file__), "data", "test")
file_list = glob.glob(os.path.join(rootdir, "*.[!p]*"))
