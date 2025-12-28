#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TSV to HTML Renderer Using Jinja2
=================================

This script renders TSV into an HTML document
using a specified Jinja2 HTML template.

Usage:
    >>> python3 html_render.py <input_tsv_file> <template_html_file>

Arguments:
    <input_tsv_file>      Path to the TSV file to render
    <template_html_file>  Path to the HTML template

Output:
    A rendered HTML file will be saved in the same location
    as the input TSV, with the same name but `.html` extension.

Requirements:
    - Jinja2 must be installed (`pip3 install jinja2`)
    - The HTML template must exist and include valid Jinja2 template variables

Author: Muhammad Ashraf
"""

import os
import sys
import pandas as pd
from jinja2 import Environment, FileSystemLoader

# Validate input arguments
if len(sys.argv) != 3:
    print("Usage: python3 html_render.py <input_tsv_file> <template_html_file>")
    sys.exit(1)

# Parse input arguments
input_tsv_file = sys.argv[1]
template_html_file = sys.argv[2]

# Generate output file name by replacing .tsv with .html
output_file = input_tsv_file.replace('.tsv', '.html')

# Extract the base name of the TSV file (without extension)
file_basename = os.path.basename(input_tsv_file).split(".")[0]

# Load TSV data into a DataFrame
df = pd.read_csv(input_tsv_file, sep='\t', keep_default_na=False, quotechar='"', dtype=str)

# Set up the Jinja2 environment using the directory of the HTML template
env = Environment(
    loader=FileSystemLoader(os.path.dirname(template_html_file)),
    autoescape=True
)

# Load the template file
template = env.get_template(os.path.basename(template_html_file))

# Render the HTML using the template and the data
rendered_html = template.render( 
    num_rows=len(df),
    tables=[df.to_html(classes='data', header="true", index=False)],
    titles=df.columns.values,
    filepath=input_tsv_file,
    file_basename=file_basename
)

rendered_html = rendered_html.replace('&lt;','<').replace('&gt;','>')

# Save the rendered HTML to file
with open(output_file, 'w', encoding="utf-8") as f:
    f.write(rendered_html)

# Print completion message with output path
print(f"HTML file has been saved to {output_file}")
