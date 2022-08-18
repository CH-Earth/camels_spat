"""Test availability of required packages."""
"""Source: https://stackoverflow.com/a/45474387"""

import os,sys
from pathlib import Path
import unittest
import pkg_resources
sys.path.append(str(Path().absolute().parent))
import python_cs_functions as cs

# Get the location of the requirements file from the config file
code_path,_ = cs.read_from_config('../0_config/config.txt','code_path')
reqs_path,_ = cs.read_from_config('../0_config/config.txt','reqs_path')
reqs_file,_ = cs.read_from_config('../0_config/config.txt','reqs_file')
_REQUIREMENTS_PATH = Path(os.path.join(code_path,reqs_path,reqs_file))

class TestRequirements(unittest.TestCase):
    """Test availability of required packages."""

    def test_requirements(self):
        """Test that each required package is available."""
        requirements = pkg_resources.parse_requirements(_REQUIREMENTS_PATH.open())
        for requirement in requirements:
            requirement = str(requirement)
            with self.subTest(requirement=requirement):
                pkg_resources.require(requirement)
                