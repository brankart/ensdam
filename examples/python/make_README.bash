#!/bin/bash

cat print_doc_*.py > make_README.py

python make_README.py > README.md
