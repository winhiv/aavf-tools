"""
Copyright Government of Canada 2018

Written by: Camy Tran, National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import os
from lxml import etree
from collections import defaultdict
import PyAAVF.parser as parser
from Asi.XML.XmlAsiTransformer import XmlAsiTransformer
from Asi.Grammar.StringMutationComparator import StringMutationComparator

"""
Takes as an input an AAVF file and outputs a csv file indicating the
drug resistance levels for each drug defined in the ASI file.
"""

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
aavf_input = TEST_PATH + '/sample.aavf'
xml_input = TEST_PATH + '/sample.xml'

def determine_resistance_levels():
    """
    Parse aavf file to AAVF object.
    Create XmlAsiTransformer.
    Use transformer to evaluate drug resistance levels in the
    AAVF object mutations.
    Output drug resistance. 
    """

    reader = parser.Reader(aavf_input)
    aavf_obj = reader.read_records()
    records = list(aavf_obj)
    mutations = defaultdict(list) # parse mutations from records

    transformer = XmlAsiTransformer(True)
    fd = open(xml_input, "r")
    genes = transformer.transform(fd)
    evaluated_genes = {}
    comparator = StringMutationComparator(True)

    output_obj = defaultdict(lambda: defaultdict(list))

    output_file = open("resistance_levels.csv","w")
    output_file.write("#gene,drug class,drug,resistance level\n")

    # create mutations list
    for record in records:
        mutations[record.GENE].append("%s%s" % (record.POS, record.REF))

    for gene in genes:
        for variant in mutations:
            evaluated_genes[gene] = genes[gene].evaluate(mutations[variant],
                                                         comparator)
    
    for gene in evaluated_genes:
        drug_classes = evaluated_genes[gene].get_evaluated_drug_classes()

        for drug_class in drug_classes:
            for drug in drug_class.get_evaluated_drugs():
                output_obj[gene][drug_class.get_drug_class()].append(drug)

    for gene in output_obj:
        for drug_class in output_obj[gene]:
            for drug in output_obj[gene][drug_class]:
                for condition in drug.get_evaluated_conditions():
                    definition = next(iter(condition.get_definitions()))
                    output_string = ("%s,%s,%s,%s\n" % (gene, drug_class.name,
                                                       drug.get_drug().name,
                                                       definition.get_text()))
                    print(output_string)
                    output_file.write(output_string)
