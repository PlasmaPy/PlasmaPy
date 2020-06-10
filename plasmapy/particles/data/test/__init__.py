import glob
import os

data_dir = os.path.dirname(__file__)
data_files = glob.glob(os.path.join(data_dir, "*.[!p]*"))
