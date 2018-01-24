#! bin/bash

message=`date`

git add *.sh 0_Scripts/*.* 0_Scripts/COG_data/*
git commit -m "Update $date"
git push origin master
