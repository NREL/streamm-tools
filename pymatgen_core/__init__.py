from __future__ import unicode_literals

__all__ = ['core','util']

import os
import warnings
import ruamel.yaml as yaml

__author__ = "Pymatgen Development Team"
__email__ ="pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ ="shyuep@gmail.com"
__version__ = "2017.8.16"

SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")

def _load_pmg_settings():
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except IOError:
        # If there are any errors, default to using environment variables
        # if present.
        d = {}
        for k, v in os.environ.items():
            if k.startswith("PMG_"):
                d[k] = v
            elif k in ["VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"]:
                d["PMG_" + k] = v
    clean_d = {}
    for k, v in d.items():
        if not k.startswith("PMG_"):
            warnings.warn('With effect from pmg 5.0, all pymatgen settings are'
                          ' prefixed with a "PMG_". E.g., "PMG_VASP_PSP_DIR" '
                          'instead of "VASP_PSP_DIR".')
            clean_d["PMG_" + k] = v
        else:
            clean_d[k] = v
    return clean_d

SETTINGS = _load_pmg_settings()


