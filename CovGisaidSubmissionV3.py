import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import argparse
from datetime import date
import glob
import shutil
import datetime
import logging
from Covid19DB import MySQLcovid19,MySQLcovid19Selector


"""
 Eric Fournier 2020-12-26
 TODO
   - faire une fonction checkup pour verifier les id manquant dans le metadata


"""


logging.basicConfig(level=logging.INFO)
pd.set_option('display.max_columns', 50)

global tosubmit_rec_list
tosubmit_rec_list = []

global tosubmit_rec_dict
tosubmit_rec_dict = {}

global covid_beluga_dir
covid_beluga_dir = "/data/PROJETS/COVID-19_Beluga/"

global covid_beluga_seq
covid_beluga_seq_dir = os.path.join(covid_beluga_dir,"Consensus/Fasta")

global covid_beluga_metadata_dir
covid_beluga_metadata_dir = os.path.join(covid_beluga_dir,"Metadata")

global published_path
published_path = "/data/PROJETS/COVID-19_Beluga/Gisaid/FinalPublished/"

global tosubmit_path
tosubmit_path = "/data/PROJETS/COVID-19_Beluga/Gisaid/ToSubmit/"

global meta_out

global mnt_partage_covid19_path
mnt_partage_covid19_path = "/mnt/Covid19/"

global partage_covid19_path
partage_covid19_path = "//swsfi52p/partage/LSPQ_Partage/DiagnosticMoleculaire/PROJETS/Covid19/"

today = date.today()
today = today.strftime("%Y%m%d")

fasta_prefix = 'hCoV-19/'

published_seqid = []

parser = argparse.ArgumentParser(description='Create files for GISAID submission')
parser.add_argument('--debug',help="run in debug mode",action='store_true')
parser.add_argument('--seq','-s',required=True,help="Fasta name having sequences to submit")
parser.add_argument('--user','-u',required=True,help="User name",choices=['foueri01','morsan01'])

args = parser.parse_args()

user = args.user
seq_in = args.seq
_debug = args.debug

if _debug:
    debug_base = "/data/Applications/GitScript/GisaidSubmission/TEST/UNPUBLISHED/"
    seq_in = os.path.join(debug_base,seq_in)
    tosubmit_path = os.path.join(debug_base,'ToSubmit')
else:
    seq_in = os.path.join(covid_beluga_seq_dir,seq_in)

if (not os.path.exists(seq_in)):
    logging.error("Input sequences or metadata file is missing")
    exit(1)

for published_fasta in glob.glob(published_path + "release1/*/all_sequences.fasta"):
    for rec in SeqIO.parse(published_fasta,'fasta'):
        rec_id = re.search(r'^hCoV-19/Canada/Qc-(\S+)/\d{4}$',rec.id).group(1)
        published_seqid.append(rec_id)

today_tosubmit_path = os.path.join(tosubmit_path,today)
meta_out = os.path.join(today_tosubmit_path,today + "_ncov19_metadata.xls")
seq_out = os.path.join(today_tosubmit_path,"all_sequences.fasta")

if not _debug:
    today_mnt_partage_covid19_path = os.path.join(mnt_partage_covid19_path,"SequencingPipeline/GISAID/TO_PUBLISH/",today)
else:
    today_mnt_partage_covid19_path = os.path.join(mnt_partage_covid19_path,"temp/",today)

try:
    os.mkdir(today_tosubmit_path)
except:
    print(today_tosubmit_path + " already exist")

def CheckMetadata(ids_in_fasta,ids_in_metadata):
    ids_in_fasta = set(ids_in_fasta)
    ids_in_metadata = set(ids_in_metadata)
    #print(ids_in_fasta)
    #print(ids_in_metadata)
    missing = ids_in_fasta - ids_in_metadata
    if len(missing) > 0:
        logging.warning("The following ids are missing from metadata : " + str(missing))
        return(True)
    return(False)
    
def UmountPartageCovid19():
    logging.info("Try to unmount Partage Covid19")
    os.system("sudo umount " + mnt_partage_covid19_path)

def MountPartageCovid19():
    logging.info("Try to mount Partage Covid19")
    UmountPartageCovid19()
    os.system("sudo mount -t cifs -o username=" + user + ",vers=3.0 " + partage_covid19_path + " " + mnt_partage_covid19_path)

UmountPartageCovid19()
MountPartageCovid19()

def from_dob_to_age(born):
    today = datetime.date.today()
    return today.year - born.year - ((today.month, today.day) < (born.month, born.day))

def GetMetadataDfFromCovBank(id_list):
    MySQLcovid19.SetConnection()
    pd_df = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),id_list)
    return(pd_df)

def BuildRecordsDictToSubmit():
    for rec in SeqIO.parse(seq_in,'fasta'):
        try:
            parsed_header = re.search(r'(Canada/Qc-)(\S+)/(\d{4}) seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+ ',rec.description)
            req_number = parsed_header.group(2)
        
            rec.id = fasta_prefix + parsed_header.group(1) + req_number + "/" + parsed_header.group(3)
            rec.description = ""
            seq_method = parsed_header.group(4)
            tosubmit_rec_list.append(rec)
            tosubmit_rec_dict[req_number] = {'gisaid_id':rec.id,'method':seq_method}
        except Exception as e:
            print(e)
            print('Unable to parse ', str(rec.id))

def GetVirusName(req_number):
    return(tosubmit_rec_dict[req_number]['gisaid_id'])

def GetSequencingMethod(req_number):
    return(tosubmit_rec_dict[req_number]['method'])

BuildRecordsDictToSubmit()
metadata_df = GetMetadataDfFromCovBank(tosubmit_rec_dict.keys())

if CheckMetadata(tosubmit_rec_dict.keys(),list(metadata_df['covv_virus_name'])):
    UmountPartageCovid19()
    exit(1)

metadata_df['covv_virus_name'] = metadata_df['covv_virus_name'].apply(GetVirusName)
metadata_df.insert(loc=11,column='covv_patient_age',value=metadata_df['DTNAISS'].apply(lambda x: from_dob_to_age(x)))
metadata_df.insert(loc=18,column='covv_seq_technology',value=metadata_df['covv_subm_sample_id'].apply(GetSequencingMethod))
del metadata_df['DTNAISS']

added_header = pd.DataFrame({'submitter':['Submitter'],'fn':['FASTA filename'],'covv_virus_name':['Virus name'],'covv_type':['Type'],'covv_passage':['Passage details/history'],'covv_collection_date':['Collection date'],'covv_location':['Location'],'covv_add_location':['Additionnal location information'],'covv_host':['Host'],'covv_add_host_info':['Additional host info'], 'covv_gender':['Gender'],'covv_patient_age':['Patient age'],'covv_patient_status':['Patient status'],'covv_specimen':['Specimen source'],'covv_outbreak':['Outbreak'],'covv_last_vaccinated':['Last vaccinated'],'covv_treatment':['Treatment'],'covv_seq_technology':['Sequencing technology'],'covv_assembly_method':['Assembly method'],'covv_coverage':['Coverage'],'covv_orig_lab':['Originating lab'],'covv_orig_lab_addr':['Address'],'covv_provider_sample_id':['Sample ID given by the sample provider'],'covv_subm_lab':['Submitting lab'],'covv_subm_lab_addr':['Address'],'covv_subm_sample_id':['Sample ID given by the submitting laboratory'],'covv_authors':['Authors']})

metadata_df = pd.concat([added_header,metadata_df])
metadata_df.to_excel(meta_out,index=False,sheet_name='Submission')
SeqIO.write(tosubmit_rec_list,seq_out,'fasta')
#print(metadata_df)

try:
    os.system("sudo mkdir {0}".format(today_mnt_partage_covid19_path))
except:
    pass
finally:
    os.system("sudo cp {0} {1} {2}".format(seq_out,meta_out,today_mnt_partage_covid19_path))

UmountPartageCovid19()

exit(0)
