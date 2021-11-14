#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "./init.sh \"Project Name\" \"CopyRight Name\" \"Author Name\" \"Link to Github Tree\""
    echo "Example:"
    echo "./init.sh \"NUBASE IAEA CRP\" \"2011, Phong\" \"Phong\" \"https://github.com/vihophong/reaclib-tools/tree/master/docs/\""
    exit 1
fi
#echo "eee"
sed -i "/project = u/c\project = u'$1'" conf.py
sed -i "/copyright = u/c\copyright = u'$2'" conf.py
sed -i "/author = u/c\author = u'$3'" conf.py
sed -i "/github_doc_root = /c\github_doc_root = '$4'" conf.py

