#!/bin/bash
#
# PY Barriat, October 2022
#
# Download and install marp (MarkDown slides extension) from here:
# https://github.com/marp-team/marp-cli/releases
#

marp --allow-local-files --theme ./assets/tum.css slides.md -o slides.pdf
marp --template bespoke --bespoke.progress --allow-local-files --theme ./assets/tum.css slides.md -o fortran.html
