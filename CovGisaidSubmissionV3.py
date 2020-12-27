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


"""
 Eric Fournier 2020-12-26
 TODO
   - utiliser le metadata au lieu de se connected sur la bd
   - faire le mount dans ce script au lieu d utiliser le sh
"""


logging.basicConfig(level=logging.INFO)
pd.set_option('display.max_columns', 50)

global qc_id_list
qc_id_list = []

global gisaid_info
gisaid_info = {}

global rec_description
rec_description = {}

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

today = date.today()
today = today.strftime("%Y%m%d")

fasta_prefix = 'hCoV-19/'

published_seqid = []


parser = argparse.ArgumentParser(description='Create files for GISAID submission')
parser.add_argument('--debug',help="run in debug mode",action='store_true')
parser.add_argument('--seq','-s',required=True,help="Fasta name having sequences to submit")

args = parser.parse_args()

seq_in = args.seq
meta_in = re.sub('^sequences','metadata',seq_in)
meta_in = re.sub('fasta$','tsv',meta_in)

_debug = args.debug

if _debug:
    debug_base = "/data/Applications/GitScript/GisaidSubmission/TEST/UNPUBLISHED/"
    seq_in = os.path.join(debug_base,seq_in)
    meta_in = os.path.join(debug_base,'metadata',meta_in)
    tosubmit_path = os.path.join(debug_base,'ToSubmit')
else:
    seq_in = os.path.join(covid_beluga_seq_dir,seq_in)
    meta_in = os.path.join(covid_beluga_metadata_dir,meta_in)

#print("SEQIN ",seq_in)
#print("METADATA ",meta_in)
#print("TOSUBMIT_PATH ",tosubmit_path)

if (not os.path.exists(seq_in)) or (not os.path.exists(meta_in)):
    logging.error("Input sequences or metadata file is missing")
    exit(1)


for published_fasta in glob.glob(published_path + "release1/*/all_sequences.fasta"):
    for rec in SeqIO.parse(published_fasta,'fasta'):
        rec_id = re.search(r'^hCoV-19/Canada/Qc-(\S+)/\d{4}$',rec.id).group(1)
        published_seqid.append(rec_id)


today_tosubmit_path = os.path.join(tosubmit_path,today)
#print("today_tosubmit_path ",today_tosubmit_path)

fasta_cat = os.path.join(today_tosubmit_path,"all_sequences.fasta")

try:
    os.mkdir(today_tosubmit_path)
except:
    print(today_tosubmit_path + " already exist")
exit(1)
#************************************************************************************
if not _debug:
    pass
    #parser.add_argument('--unpublished-path','-u', required=True, help="path to final unpublished sequence")
    #parser.add_argument('--submitted-path','-s', required=True, help="path to final sumitted sequence")

    #args = parser.parse_args()

    #unpublished_path = os.path.join(args.unpublished_path)
    #today_submitted_seq_path = os.path.join(args.submitted_path,today)
    #metadata_name = today + "_ncov19_metadata.xls"

    #lspq_miseq_dir = os.path.join("/mnt/Partage/LSPQ_Partage/DiagnosticMoleculaire/PROJETS/Covid19/SequencingPipeline/GISAID/TO_PUBLISH")
else:
    pass
    #unpublished_path = "/data/Applications/GitScript/GisaidSubmission/TEST/UNPUBLISHED/"
    #today_submitted_seq_path = "/data/Applications/GitScript/GisaidSubmission/TEST/FINAL_SUBMITTED/"
    #metadata_name = today + "_ncov19_metadata.xls"

    #lspq_miseq_dir = "/data/Applications/GitScript/GisaidSubmission/TEST/TO_PUBLISH/"





try:
    os.mkdir(today_submitted_seq_path)
except:
    print(today_submitted_seq_path + " already exist")


def GetSequencingMethod(samples_id):
    method_list = []
    for sample in samples_id:
         try:
             method = re.search(r'seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+',str(gisaid_info[sample]['method'])).group(1)
         except:
             method = 'unknown'
         
         method_list.append(method)
    
    return(method_list)

def ConcatSeq(files_in,file_out):
    rec_list = []
    for fasta in glob.glob(files_in):
        my_rec = SeqIO.read(fasta,'fasta')
        my_rec.id = fasta_prefix + my_rec.id
        desc = re.split('\s+',my_rec.description)[1]
        my_rec.description = ''
        rec_list.append(my_rec)
        rec_description[my_rec.id] = desc
    SeqIO.write(rec_list,file_out,'fasta')

def MakeSeqIdList():

    for rec in SeqIO.parse(fasta_cat,'fasta'):
        try:
            #TODO enlever le prefix L du regex
            qc_id = re.search(fasta_prefix + r'Canada/Qc-(L\S+)/\d{4}',rec.id).group(1)
            if qc_id not in published_seqid:
                qc_id_list.append(qc_id)
            
            desc_dict = {}
            desc_dict['gisaid_id'] = rec.id
            desc_dict['method'] = rec_description[rec.id]
            
            gisaid_info[qc_id] = desc_dict

        except Exception as e:
            print(e)
            print('Unable to parse ', str(rec.id))


def from_dob_to_age(born):
    today = datetime.date.today()
    return today.year - born.year - ((today.month, today.day) < (born.month, born.day))

def GetQcDataframeFromDSPdb(id_list):
    MySQLcovid19.SetConnection()
    pd_df = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),id_list)
    return(pd_df)

ConcatSeq(unpublished_path + "*.fasta",fasta_cat)

MakeSeqIdList()

'''
def Test(qc_list):
    MySQLcovid19.SetConnection()
    for spec_id in qc_list:
        pd_df = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),[spec_id])
        #print(spec_id)
Test(qc_id_list)
'''

df_qc_from_dsp_db = GetQcDataframeFromDSPdb(qc_id_list)
df_qc_from_dsp_db['covv_virus_name'] = 'hCoV-19/Canada/Qc-' +  df_qc_from_dsp_db['covv_virus_name'] + "/2020"
df_qc_from_dsp_db.insert(loc=11,column='covv_patient_age',value=df_qc_from_dsp_db['DTNAISSINFO'].apply(lambda x: from_dob_to_age(x)))
df_qc_from_dsp_db.insert(loc=18,column='covv_seq_technology',value= GetSequencingMethod(df_qc_from_dsp_db['covv_subm_sample_id']),allow_duplicates=True)

del df_qc_from_dsp_db['DTNAISSINFO']

added_header = pd.DataFrame({'submitter':['Submitter'],'fn':['FASTA filename'],'covv_virus_name':['Virus name'],'covv_type':['Type'],'covv_passage':['Passage details/history'],'covv_collection_date':['Collection date'],'covv_location':['Location'],'covv_add_location':['Additionnal location information'],'covv_host':['Host'],'covv_add_host_info':['Additional host info'], 'covv_gender':['Gender'],'covv_patient_age':['Patient age'],'covv_patient_status':['Patient status'],'covv_specimen':['Specimen source'],'covv_outbreak':['Outbreak'],'covv_last_vaccinated':['Last vaccinated'],'covv_treatment':['Treatment'],'covv_seq_technology':['Sequencing technology'],'covv_assembly_method':['Assembly method'],'covv_coverage':['Coverage'],'covv_orig_lab':['Originating lab'],'covv_orig_lab_addr':['Address'],'covv_provider_sample_id':['Sample ID given by the sample provider'],'covv_subm_lab':['Submitting lab'],'covv_subm_lab_addr':['Address'],'covv_subm_sample_id':['Sample ID given by the submitting laboratory'],'covv_authors':['Authors']})


df_qc_from_dsp_db = pd.concat([added_header,df_qc_from_dsp_db])
print(df_qc_from_dsp_db.shape)
metadata_out = os.path.join(today_submitted_seq_path,metadata_name)

df_qc_from_dsp_db.to_excel(metadata_out,index=False,sheet_name='Submissions')

fasta_in_list = glob.glob(unpublished_path + "*.fasta")

for fasta in fasta_in_list:
    shutil.move(fasta,today_submitted_seq_path)

os.system("sudo mkdir " + os.path.join(lspq_miseq_dir,today))
os.system("sudo cp " + fasta_cat + " " + metadata_out + " " + os.path.join(lspq_miseq_dir,today))

exit(0)
