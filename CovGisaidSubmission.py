import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


"""
Eric Fournier 2020-05-14


"""



pd.set_option('display.max_columns', 50)

global sgil_id_list
sgil_id_list = []

global gisaid_info
gisaid_info = {}

base_dir = "/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/FinalSubmitted_20200514/"
fasta_cat = os.path.join(base_dir,"all_sequences.fasta")
metadata_out = os.path.join(base_dir,"20200514_ncov19_metadata.xls")
submitting_lab_addr = "20045, chemin Sainte-Marie, Sainte-Anne-de-Bellevue, QC, Canada"
authors = "Sandrine Moreira, Ioannis Ragoussis, Guillaume Bourque, Jesse Shapiro, Mark Lathrop and Michel Roger on behalf of the CoVSeQ research group (http://covseq.ca/researchgroup)"

def ConcatSeq(files_in,file_out):
    os.system("cat " + files_in + " > " + file_out )

def MakeSeqIdList():

    for rec in SeqIO.parse(fasta_cat,'fasta'):
        try:
            sgil_id = re.search(r'^Canada/Qc-(L\S+)/\d{4}',rec.id).group(1)
            sgil_id_list.append(sgil_id)
            
            desc_dict = {}
            desc_dict['gisaid_id'] = rec.id
            desc_dict['method'] = re.split('\s+',rec.description)[1]
            
            gisaid_info[sgil_id] = desc_dict

        except:
            print('Unable to parse ', str(rec.id))


def GetSequencingMethod(samples_id):
    method_list = []
    for sample in samples_id:
         try:
             method = re.search(r'seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+',str(gisaid_info[sample]['method'])).group(1)
         except:
             method = 'unknown'
         
         method_list.append(method)
    
    return(method_list)

#ConcatSeq(base_dir + "L*.fasta",fasta_cat)

MakeSeqIdList()

in_metadata = os.path.join(base_dir,"sgil_extract.tsv")
df_in = pd.read_csv(in_metadata,delimiter="\t",index_col=False,encoding="UTF-8")
sub_df_in = df_in[df_in['NO_LSPQ'].isin(sgil_id_list)]
sub_df_in.insert(loc=0,column='ID_GISAID',value= "Canada/Qc-" + sub_df_in['NO_LSPQ'] + "/2020" ,allow_duplicates=False ) 
sub_df_in.insert(loc=1,column='LOCATION_GISAID',value= "North-America / Canada / " + sub_df_in['RSS_PATIENT'] ,allow_duplicates=False ) 
sub_df_in.insert(loc=2,column='TRAVEL_GISAID',value= sub_df_in['VOYAGE_PAYS_1'] ,allow_duplicates=False ) 

gisaid_metadata = sub_df_in.loc[:,('NO_LSPQ','ID_GISAID','DATE_PRELEV','AGE','SEX','LOCATION_GISAID','TRAVEL_GISAID','RES_LAB','CH','ADDR_CH')]

gisaid_metadata.columns = ['NO_LSPQ temp','Virus name temp','Collection date temp','Patient age temp','Gender temp','Location temp','Additional location information temp','Additional host information temp','Originating lab temp','Originating lab address temp']

gisaid_metadata.insert(loc=0,column='NO_LSPQ',value=gisaid_metadata['NO_LSPQ temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=1,column='submitter',value="SandrineMoreiraLSPQ",allow_duplicates=True)
gisaid_metadata.insert(loc=2,column='fn',value=os.path.basename(fasta_cat),allow_duplicates=True)
gisaid_metadata.insert(loc=3,column='covv_virus_name',value=gisaid_metadata['Virus name temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=4,column='covv_type',value='betacoronavirus',allow_duplicates=True)
gisaid_metadata.insert(loc=5,column='covv_passage',value='Original',allow_duplicates=True)
gisaid_metadata.insert(loc=6,column='covv_collection_date',value=gisaid_metadata['Collection date temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=7,column='covv_location',value=gisaid_metadata['Location temp'],allow_duplicates=True)
#gisaid_metadata.insert(loc=8,column='covv_add_location',value=gisaid_metadata['Additional location information temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=8,column='covv_add_location',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=9,column='covv_host',value='Human',allow_duplicates=True)
#gisaid_metadata.insert(loc=10,column='covv_add_host_info',value=gisaid_metadata['Additional host information temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=10,column='covv_add_host_info',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=11,column='covv_gender',value=gisaid_metadata['Gender temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=12,column='covv_patient_age',value=gisaid_metadata['Patient age temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=13,column='covv_patient_status',value='unknown',allow_duplicates=True)
gisaid_metadata.insert(loc=14,column='covv_specimen',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=15,column='covv_outbreak',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=16,column='covv_last_vaccinated',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=17,column='covv_treatment',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=18,column='covv_seq_technology',value= GetSequencingMethod(gisaid_metadata['NO_LSPQ']),allow_duplicates=True)
gisaid_metadata.insert(loc=19,column='covv_assembly_method',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=20,column='covv_coverage',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=21,column='covv_orig_lab',value=gisaid_metadata['Originating lab temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=22,column='covv_orig_lab_addr',value=gisaid_metadata['Originating lab address temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=23,column='covv_provider_sample_id',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=24,column='covv_subm_lab',value="Laboratoire de santé publique du Québec",allow_duplicates=True)
gisaid_metadata.insert(loc=25,column='covv_subm_lab_addr',value=submitting_lab_addr,allow_duplicates=True)
gisaid_metadata.insert(loc=26,column='covv_subm_sample_id',value=gisaid_metadata['NO_LSPQ temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=27,column='covv_authors',value=authors,allow_duplicates=True)


del gisaid_metadata['NO_LSPQ']
del gisaid_metadata['NO_LSPQ temp']
del gisaid_metadata['Virus name temp']
del gisaid_metadata['Collection date temp']
del gisaid_metadata['Patient age temp']
del gisaid_metadata['Gender temp']
del gisaid_metadata['Location temp']
del gisaid_metadata['Additional location information temp']
del gisaid_metadata['Additional host information temp']
del gisaid_metadata['Originating lab temp']
del gisaid_metadata['Originating lab address temp']
#print(gisaid_metadata)

added_header = {'submitter':'Submitter','fn':'FASTA filename','covv_virus_name':'Virus name','covv_type':'Type','covv_passage':'Passage details/history','covv_collection_date':'Collection date','covv_location':'Location','covv_add_location':'Additionnal location information','covv_host':'Host','covv_add_host_info':'Additional host info', 'covv_gender':'Gender','covv_patient_age':'Patient age','covv_patient_status':'Patient status','covv_specimen':'Specimen source','covv_outbreak':'Outbreak','covv_last_vaccinated':'Last vaccinated','covv_treatment':'Treatment','covv_seq_technology':'Sequencing technology','covv_assembly_method':'Assembly method','covv_coverage':'Coverage','covv_orig_lab':'Originating lab','covv_orig_lab_addr':'Address','covv_provider_sample_id':'Sample ID given bt the sample provider','covv_subm_lab':'Submitting lab','covv_subm_lab_addr':'Address','covv_subm_sample_id':'Sample ID given by the submitting laboratory','covv_authors':'Authors'}

gisaid_metadata.iloc[0] = added_header

gisaid_metadata.to_excel(metadata_out,index=False,sheet_name='Submissions')
