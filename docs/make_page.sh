#!/bin/bash

make clean
make html
rsync -av _build/html/* ../html/
