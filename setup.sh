#!/bin/bash

#0. Install singularity
which singularity && singularity pull docker://benpasto/smrnaseq

#1. Install nextflow
which nextflow || curl -s https://get.nextflow.io | bash

if [[ ":$PATH:" == *":$PWD:"* ]]; then
	echo "in path"
else
	echo "not"
	echo "export PATH=\$PATH:\"$PWD\"" >> ~/.bashrc
	source ~/.bashrc
fi
