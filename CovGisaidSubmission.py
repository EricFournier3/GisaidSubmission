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


pd.set_option('display.max_columns', 20)

base_dir = "/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/FinalSubmitted_20200514/"
fasta_cat = os.path.join(base_dir,"all_sequences.fasta")

def ConcatSeq(files_in,file_out):
    os.system("cat " + files_in + " > " + file_out )

def MakeSeqIdList():
    sgil_id_list = []
    gisaid_info = {}
    #Canada/Qc-L00240556/2020
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

    print(str(gisaid_info))
    print(str(sgil_id_list))

MakeSeqIdList()

exit(0)

ConcatSeq(base_dir + "*.fasta",fasta_cat)



in_metadata = os.path.join(base_dir,"20200508_SelectionPatients.txt")

df_in = pd.read_csv(in_metadata,delimiter="\t",index_col=False,encoding="ISO-8859-1")
df_in["ID_GISAID"] = "Canada/Qc-" + df_in["sample"] + "/2020"
df_in["LOCATION_GISAID"] = "North-America / Canada / Quebec / " + df_in[" RSS_PATIENT"]
df_in["TRAVEL_GISAID"] = df_in[" VOYAGE_PAYS_1"]

gisaid_metadata = df_in.loc[:,('ID_GISAID',' DATE_PRELEVEMENT',' AGE_ANNEE',' SEXE','LOCATION_GISAID','TRAVEL_GISAID',' RESULTAT_LABORATOIRE',' CH')]

#print(gisaid_metadata)

gisaid_metadata.columns = ['Virus name temp','Collection date temp','Patient age temp','Gender temp','Location temp','Additional location information temp','Additional host information temp','Originating lab temp']

gisaid_metadata.insert(loc=0,column='Submitter',value="SandrineMoreiraLSPQ",allow_duplicates=True)
gisaid_metadata.insert(loc=1,column='FASTA filename',value=os.path.basename(fasta_cat),allow_duplicates=True)
gisaid_metadata.insert(loc=2,column='Virus name',value=gisaid_metadata['Virus name temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=3,column='Type',value='betacoronavirus',allow_duplicates=True)
gisaid_metadata.insert(loc=4,column='Passage details/history',value='Original',allow_duplicates=True)
gisaid_metadata.insert(loc=5,column='Collection date',value=gisaid_metadata['Collection date temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=6,column='Location',value=gisaid_metadata['Location temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=7,column='Additional location information',value=gisaid_metadata['Additional location information temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=8,column='Host',value='Human',allow_duplicates=True)
gisaid_metadata.insert(loc=9,column='Additional host information',value=gisaid_metadata['Additional host information temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=10,column='Gender',value=gisaid_metadata['Gender temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=11,column='Patient age',value=gisaid_metadata['Patient age temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=12,column='Patient status',value='unknown',allow_duplicates=True)
gisaid_metadata.insert(loc=13,column='Specimen source',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=14,column='Outbreak',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=15,column='Last vaccinated',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=16,column='Treatment',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=17,column='Sequencing technology',value='MGI',allow_duplicates=True)
gisaid_metadata.insert(loc=18,column='Assembly method',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=19,column='Coverage',value=' ',allow_duplicates=True)
gisaid_metadata.insert(loc=20,column='Originating lab',value=gisaid_metadata['Originating lab temp'],allow_duplicates=True)
gisaid_metadata.insert(loc=21,column='Originating lab address',value="a venir",allow_duplicates=True)
gisaid_metadata.insert(loc=22,column='Sample ID given by the sample provider',value=" ",allow_duplicates=True)
gisaid_metadata.insert(loc=23,column='Submitting lab',value="Laboratoire de santé publique du Québec",allow_duplicates=True)
gisaid_metadata.insert(loc=24,column='Submitting lab address',value="a venir",allow_duplicates=True)
gisaid_metadata.insert(loc=25,column='Sample ID given by the submitting laboratory',value=" ",allow_duplicates=True)
gisaid_metadata.insert(loc=26,column='Authors',value="a venir",allow_duplicates=True)

print(gisaid_metadata)


