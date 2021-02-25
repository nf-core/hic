#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# Add additional regexes for new tools in process get_software_versions
regexes = {
    'nf-core/hic': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'Bowtie2': ['v_bowtie2.txt', r"bowtie2-align-s version (\S+)"],
    'Python': ['v_python.txt', r"Python (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['nf-core/hic'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie2'] = '<span style="color:#999999;\">N/A</span>'
results['Python'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/hic Software Versions'
section_href: 'https://github.com/nf-core/hic'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
