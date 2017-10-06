import os


HOME_DIR = os.getcwd()
RELATIVE_TEST_DIR = os.path.join(os.path.dirname(__file__))
TEST_DIR = os.path.join(HOME_DIR,RELATIVE_TEST_DIR)

os.chdir(TEST_DIR)
