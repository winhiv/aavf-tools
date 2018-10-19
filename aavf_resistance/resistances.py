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
from PyAAVF import parser
from Asi.XML.XmlAsiTransformer import XmlAsiTransformer
from Asi.Grammar.StringMutationComparator import StringMutationComparator

# Takes as an input an AAVF file and outputs a csv file indicating the
# drug resistance levels for each drug defined in the input ASI file.

# default input and output values
DEFAULT_PATH = os.path.dirname(os.path.abspath(__file__))


@click.command()
@click.option('-a', '--aavf_input',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-x', '--xml_input',
              type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option('-o', '--output', type=click.File('w'))
def determine_resistance_levels(aavf_input, xml_input, output):
    """
    An example tool that queries an aavf file against
    Stanford HIVdb to determine drug resistance levels
    for drugs defined in an ASI file.

    By default, this example tool will run on the sample files
    in aavf_resistance/data, and output to
    aavf_resistance/output/resistance_levels.csv
    """

    aavf_file = DEFAULT_PATH + '/data/sample.aavf'
    xml_file = DEFAULT_PATH + '/data/sample.xml'
    output_path = DEFAULT_PATH + '/output/resistance_levels.csv'
    # default input and output values

    if aavf_input:
        aavf_file = aavf_input

    if xml_input:
        xml_file = xml_input

    if output:
        output_path = output

    output_string = output_resistance_levels(aavf_file, xml_file, output_path)

    click.echo(output_string)


def output_resistance_levels(aavf_file, xml_file, output_path):
    """
    Parse aavf file to AAVF object.
    Create XmlAsiTransformer.
    Use transformer to evaluate drug resistance levels in the
    AAVF object mutations.
    Output drug resistance.
    """

    records = list(parser.Reader(aavf_file).read_records())
    mutations = defaultdict(list)  # parse mutations from records
    genes = XmlAsiTransformer(True).transform(open(xml_file, "r"))

    output_file = output_path

    # handles the case where the output_path was passed from somewhere
    # other than to click, from the command line,
    # e.g. from a test file
    if isinstance(output_path, str):
        output_file = open(output_path, "w+")

    output_string = "#gene,drug class,drug,resistance level"

    # create mutations list
    for record in records:
        if record.ALT[0] != "*":
            mutations[record.GENE].append("%s%s" % (record.POS, record.ALT[0]))

    for gene in mutations:
        evaluated_gene = genes[gene].evaluate(mutations[gene],
                                              StringMutationComparator(True))

        for drug_class in evaluated_gene.get_evaluated_drug_classes():
            for drug in drug_class.get_evaluated_drugs():
                for condition in drug.get_evaluated_conditions():
                    definition = next(iter(condition.get_definitions()))

                    output_string += ("\n%s,%s,%s,%s" %
                                      (gene, drug_class.get_drug_class().name,
                                       drug.get_drug().name,
                                       definition.get_text()))

    output_file.write(output_string)
    output_file.close()

    return output_string
