#!/bin/bash
bedtools getfasta -fi $1 -bed $2 -fo $3 -name
echo "Finished Promotor Sequence Extraction"
