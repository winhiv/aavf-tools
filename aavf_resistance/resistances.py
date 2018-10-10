"""
Copyright Government of Canada 2018

Written by: Camy Tran, National Microbiology Laboratory, Public Health Agency
of Canada

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
from collections import defaultdict
import click
import PyAAVF.parser as parser
from Asi.XML.XmlAsiTransformer import XmlAsiTransformer
from Asi.Grammar.StringMutationComparator import StringMutationComparator


"""
Takes as an input an AAVF file and outputs a csv file indicating the
drug resistance levels for each drug defined in the input ASI file.
"""

# default input and output values
DEFAULT_PATH = os.path.dirname(os.path.abspath(__file__))


@click.command()
@click.option('-a', '--aavf_input',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-x', '--xml_input',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-o', '--output', type=click.File('w'))
def determine_resistances(aavf_input, xml_input, output):
    """
    An example tool that queries an aavf file against
    Stanford HIVdb to determine drug resistance levels
    for drugs defined in an ASI file.

    By default, this example tool will run on the sample files
    in aavf_resistance/data, and output to
    aavf_resistance/output/resistance_levels.csv

    Parse aavf file to AAVF object.
    Create XmlAsiTransformer.
    Use transformer to evaluate drug resistance levels in the
    AAVF object mutations.
    Output drug resistance.
    """

    aavf_file = DEFAULT_PATH + '/data/sample.aavf'
    xml_file = DEFAULT_PATH + '/data/sample.xml'
    output_path = DEFAULT_PATH + '/output/resistance_levels.csv'
    # default input and output values

    if aavf_input:
        aavf_file = aavf_input

    if xml_input:
        xml_file = xml_input

    reader = parser.Reader(aavf_file)
    records = list(reader.read_records())
    mutations = defaultdict(list)  # parse mutations from records

    transformer = XmlAsiTransformer(True)
    genes = transformer.transform(open(xml_file, "r"))
    comparator = StringMutationComparator(True)

    if output:
        output_file = output
    else:
        output_file = open(output_path, "w+")

    output_file.write("#gene,drug class,drug,resistance level\n")

    # create mutations list
    for record in records:
        if record.ALT[0] != "*":
            mutations[record.GENE].append("%s%s" % (record.POS, record.ALT[0]))

    for gene in mutations:
        evaluated_gene = genes[gene].evaluate(mutations[gene], comparator)
        drug_classes = evaluated_gene.get_evaluated_drug_classes()

        for drug_class in drug_classes:
            for drug in drug_class.get_evaluated_drugs():
                for condition in drug.get_evaluated_conditions():
                    definition = next(iter(condition.get_definitions()))

                    output_string = ("%s,%s,%s,%s\n" %
                                     (gene, drug_class.get_drug_class().name,
                                      drug.get_drug().name,
                                      definition.get_text()))

                    click.echo(output_string)
                    output_file.write(output_string)

    output_file.close()
