#!/bin/bash

echo "TESTING"

unpublished_path="/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/TEST/FinalUnpublished/"
submitted_path="/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/TEST/FinalSubmitted/"


<<COM
nb_mounted_rep=$(ls /mnt/Partage/ | wc -l)

if [ $nb_mounted_rep -gt 0 ]
        then
        echo -e "${green_message}INFO: " "/mnt/Partage/ already mounted"
        echo -e "${white_message}"

else
        echo -e "${yellow_message}INFO: " "/mnt/Partage/ not mounted. Try to mount now ..."
        echo -e "${white_message}"
        read pw < $PASS_FILE
        sudo mount -t cifs -o username=${ldap_user},password=$pw,vers=3.0 "//swsfi52p/partage" /mnt/Partage
fi
COM




#python /data/Applications/GitScript/GisaidSubmission/CovGisaidSubmission.py --unpublished-path  ${unpublished_path} --submitted-path ${submitted_path}
