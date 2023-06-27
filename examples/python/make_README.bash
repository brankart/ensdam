#!/bin/bash

cat print_header_doc.py > make_README.py
cat print_doc_*.py >> make_README.py

python make_README.py > README.md
