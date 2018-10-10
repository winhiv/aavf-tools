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

from setuptools import find_packages, setup

dependencies = ['click', "PyAAVF"]

setup(
    name='aavf-tools',
    version='0.0.0',
    url='https://github.com/winhiv/aavf-tools.git',
    license='Apache License, Version 2.0',
    author='Camy Tran',
    author_email='camy.tran@canada.ca',
    description='A collection of example tools for annotating/querying aavf files against stanford hivdb',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    dependency_links=['git+https://github.com/winhiv/PyAAVF#egg=PyAAVF'],
    entry_points='''
        [console_scripts]
        aavfresistance=aavf_resistance.resistances:determine_resistance_levels

    ''',
    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Environment :: Console',
        'Programming Language :: Python :: 3',
    ]
)
