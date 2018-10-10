"""
Copyright Government of Canada 2018

Written by: Camy Tran, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not us
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import os
from Asi.Grammar.StringMutationComparator import StringMutationComparator
from aavf_resistance.resistances import output_resistance_levels

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
AAVF_FILE = TEST_PATH + '/data/sample.aavf'
XML_FILE = TEST_PATH + '/data/sample.xml'

OUTPUT_FILE = TEST_PATH + '/output/sample_resistance_levels.csv'

class TestAAVFResistance:

    @classmethod
    def setup_class(cls):
        cls.strict_comparison = True
        cls.mutation_comparator = StringMutationComparator(cls.strict_comparison)

    def test_output_file(self):
        output_resistance_levels(AAVF_FILE, XML_FILE, OUTPUT_FILE)    

        assert os.path.isfile(AAVF_FILE)
        assert os.path.isfile(OUTPUT_FILE)

        if os.path.isfile(OUTPUT_FILE):
            os.remove(OUTPUT_FILE)
