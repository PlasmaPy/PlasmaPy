"""
Test`.h5` files came from
https://github.com/openPMD/openPMD-example-datasets/tree/b4f87b817629b99a048026a2724a5b265810d8be
"""
import glob
import os

data_dir = os.path.dirname(__file__)
data_files = glob.glob(os.path.join(data_dir, "*.[!p]*"))
