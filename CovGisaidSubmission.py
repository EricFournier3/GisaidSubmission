import pandas as pd
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


"""
Eric Fournier 2020-05-14


"""

#TODO concat fasta
#TODO Extraire sequencing technologie a partir du fasta
#TODO extraire du metadata seulement les specimens qui sont dans le fasta


pd.set_option('display.max_columns', 50)

global sgil_id_list
sgil_id_list = []

global gisaid_info
gisaid_info = {}

base_dir = "/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/FinalSubmitted_20200514/"
fasta_cat = os.path.join(base_dir,"all_sequences.fasta")

def ConcatSeq(files_in,file_out):
    os.system("cat " + files_in + " > " + file_out )

def MakeSeqIdList():

    for rec in SeqIO.parse(fasta_cat,'fasta'):
        try:
            sgil_id = re.search(r'^Canada/Qc-(L\S+)/\d{4}',rec.id).group(1)
            sgil_id_list.append(sgil_id)
            
            #print(rec.id, "  description ", rec.description)
            desc_dict = {}
            desc_dict['gisaid_id'] = rec.id
            desc_dict['method'] = re.split('\s+',rec.description)[1]
            
            gisaid_info[sgil_id] = desc_dict

        except:
            print('Unable to parse ', str(rec.id))



#ConcatSeq(base_dir + "L*.fasta",fasta_cat)

MakeSeqIdList()

#print(str(gisaid_info))

#in_metadata = os.path.join(base_dir,"20200508_SelectionPatients.txt")
in_metadata = os.path.join(base_dir,"sgil_extract.tsv")

df_in = pd.read_csv(in_metadata,delimiter="\t",index_col=False,encoding="UTF-8")
#df_in = pd.read_csv(in_metadata,delimiter="\t",index_col=False,encoding="iso-8859-1")


#sub_df_in = df_in[df_in['sample'].isin(sgil_id_list)]
sub_df_in = df_in[df_in['NO_LSPQ'].isin(sgil_id_list)]

#sub_df_in["ID_GISAID"] = "Canada/Qc-" + sub_df_in["sample"] + "/2020"
#sub_df_in["ID_GISAID"] = "Canada/Qc-" + sub_df_in["NO_LSPQ"] + "/2020"

sub_df_in.insert(loc=0,column='ID_GISAID',value= "Canada/Qc-" + sub_df_in['NO_LSPQ'] + "/2020" ,allow_duplicates=False ) 

#sub_df_in["LOCATION_GISAID"] = "North-America / Canada / Quebec / " + sub_df_in[" RSS_PATIENT"]
#sub_df_in["LOCATION_GISAID"] = "North-America / Canada / Quebec / " + sub_df_in["RSS_PATIENT"]
sub_df_in.insert(loc=1,column='LOCATION_GISAID',value= "North-America / Canada / " + sub_df_in['RSS_PATIENT'] ,allow_duplicates=False ) 

#sub_df_in["TRAVEL_GISAID"] = sub_df_in[" VOYAGE_PAYS_1"]
#sub_df_in["TRAVEL_GISAID"] = sub_df_in["VOYAGE_PAYS_1"]
sub_df_in.insert(loc=2,column='TRAVEL_GISAID',value= sub_df_in['VOYAGE_PAYS_1'] ,allow_duplicates=False ) 

#print(sub_df_in)

#gisaid_metadata = sub_df_in.loc[:,('sample','ID_GISAID',' DATE_PRELEVEMENT',' AGE_ANNEE',' SEXE','LOCATION_GISAID','TRAVEL_GISAID',' RESULTAT_LABORATOIRE',' CH',)]

gisaid_metadata = sub_df_in.loc[:,('NO_LSPQ','ID_GISAID','DATE_PRELEV','AGE','SEX','LOCATION_GISAID','TRAVEL_GISAID','RES_LAB','CH','ADDR_CH')]

gisaid_metadata.columns = ['NO_LSPQ temp','Virus name temp','Collection date temp','Patient age temp','Gender temp','Location temp','Additional location information temp','Additional host information temp','Originating lab temp','Originating lab address temp']


def GetSequencingMethod(samples_id):
    method_list = []
    for sample in samples_id:
         try:
             method = re.search(r'seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+',str(gisaid_info[sample]['method'])).group(1)
         except:
             method = 'unknown'
         
         method_list.append(method)
    
    return(method_list)


gisaid_metadata.insert(loc=0,column='NO_LSPQ',value=gisaid_metadata['NO_LSPQ temp'],allow_duplicates=True)
#myt = test(gisaid_metadata['NO_LSPQ'].values)
#print(str(myt))
gisaid_metadata.insert(loc=1,column='Submitter',value="SandrineMoreiraLSPQ",allow_duplicates=True)
gisaid_metadata.insert(loc=2,column='FASTA filename',value=os.path.basename(fasta_cat),allow_duplicates=True)
gisaid_metadata.insert(loc=3,column='Virus name',value=gisaid_metadata['Virus name temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=4,column='Type',value='betacoronavirus',allow_duplicates=True)
gisaid_metadata.insert(loc=5,column='Passage details/history',value='Original',allow_duplicates=True)
gisaid_metadata.insert(loc=6,column='Collection date',value=gisaid_metadata['Collection date temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=7,column='Location',value=gisaid_metadata['Location temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=8,column='Additional location information',value=gisaid_metadata['Additional location information temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=9,column='Host',value='Human',allow_duplicates=True)
gisaid_metadata.insert(loc=10,column='Additional host information',value=gisaid_metadata['Additional host information temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=11,column='Gender',value=gisaid_metadata['Gender temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=12,column='Patient age',value=gisaid_metadata['Patient age temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=13,column='Patient status',value='unknown',allow_duplicates=True)
gisaid_metadata.insert(loc=14,column='Specimen source',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=15,column='Outbreak',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=16,column='Last vaccinated',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=17,column='Treatment',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=18,column='Sequencing technology',value= GetSequencingMethod(gisaid_metadata['NO_LSPQ']),allow_duplicates=True)
#gisaid_metadata.insert(loc=18,column='Sequencing technology',value= "jjjjj",allow_duplicates=True)
gisaid_metadata.insert(loc=19,column='Assembly method',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=20,column='Coverage',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=21,column='Originating lab',value=gisaid_metadata['Originating lab temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=22,column='Originating lab address',value=gisaid_metadata['Originating lab address temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=23,column='Sample ID given by the sample provider',value=" ",allow_duplicates=True)
gisaid_metadata.insert(loc=24,column='Submitting lab',value="Laboratoire de santé publique du Québec",allow_duplicates=True)
gisaid_metadata.insert(loc=25,column='Submitting lab address',value="a venir",allow_duplicates=True)
gisaid_metadata.insert(loc=26,column='Sample ID given by the submitting laboratory',value=" ",allow_duplicates=True)
gisaid_metadata.insert(loc=27,column='Authors',value="a venir",allow_duplicates=True)


del gisaid_metadata['NO_LSPQ temp']
#del gisaid_metadata['NO_LSPQ']
del gisaid_metadata['Virus name temp']
del gisaid_metadata['Collection date temp']
del gisaid_metadata['Patient age temp']
del gisaid_metadata['Gender temp']
del gisaid_metadata['Location temp']
del gisaid_metadata['Additional location information temp']
del gisaid_metadata['Additional host information temp']
del gisaid_metadata['Originating lab temp']
del gisaid_metadata['Originating lab address temp']
print(gisaid_metadata)


