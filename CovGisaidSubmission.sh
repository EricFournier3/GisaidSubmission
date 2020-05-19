#!/bin/bash

echo "TESTING"

unpublished_path="/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/TEST/FinalUnpublished/"
submitted_path="/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/TEST/FinalSubmitted/"

slbio_user=$(whoami)
ldap_user=$(echo ${slbio_user} | cut -d '@' -f 1)
PASS_FILE="/home/${slbio_user}/pass.txt"

nb_mounted_rep=$(ls /mnt/Partage/ | wc -l)

green_message="\e[32m"
white_message="\e[39m"
red_message="\e[31m"
yellow_message="\e[33m"


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


#python /data/Applications/GitScript/GisaidSubmission/CovGisaidSubmission.py --unpublished-path  ${unpublished_path} --submitted-path ${submitted_path}
