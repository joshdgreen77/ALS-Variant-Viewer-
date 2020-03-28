#!/bin/bash
for i in $(ls -1 *.csv); do
	rscript variant-fixer.R $i 
done