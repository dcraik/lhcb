#!/bin/bash
m4 -DM4_DATAYEAR=2016 -DTEST -P Run2Jets_zjet.template.py > Run2Jets_zjet_test.py 
m4 -DM4_DATAYEAR=2015 -P Run2Jets_zjet.template.py > Run2Jets_zjet15.py 
m4 -DM4_DATAYEAR=2016 -P Run2Jets_zjet.template.py > Run2Jets_zjet16.py 
m4 -DM4_DATAYEAR=2017 -P Run2Jets_zjet.template.py > Run2Jets_zjet17.py 
m4 -DM4_DATAYEAR=2018 -P Run2Jets_zjet.template.py > Run2Jets_zjet18.py 
