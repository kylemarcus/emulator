#! /bin/bash

cd hadoop-0.20.1/
./bin/hadoop namenode -format
./bin/start-dfs.sh

