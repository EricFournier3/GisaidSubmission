# coding=utf-8
"""
Eric Fournier 2020-08-12

"""

import mysql.connector
import pandas as pd

class MySQLcovid19:
    host = 'localhost'
    user = 'root'
    password = 'lspq2019'
    database = 'TestCovid19_20200911'
    connection = None

    @classmethod
    def SetConnection(cls):
        cls.connection =  mysql.connector.connect(host=cls.host,user=cls.user,password=cls.password,database=cls.database)
  
    @classmethod
    def GetConnection(cls):
        return cls.connection

    @classmethod
    def GetCursor(cls):
        return cls.GetConnection().cursor()

    @classmethod
    def Commit(cls):
        cls.connection.commit()


class MySQLcovid19Selector:

    @classmethod
    def GetSampleDate(cls,cursor,sample_name):

        cursor.execute("select DATE_PRELEV_HOPITAL from Prelevements where GENOME_QUEBEC_REQUETE  like '%{0}%'".format(sample_name))
        #date_prelev = cursor.fetchone()
        date_prelev = cursor.fetchall()

        if len(date_prelev) == 1 :
            date_prelev = date_prelev[0][0]
            mois = Utils.GetFrenchMonth(date_prelev.strftime("%B"))
            return(date_prelev.strftime("%d {0} %Y".format(mois)))

        return("indéterminé")

    @classmethod
    def GetMetadataAsPdDataFrame(cls,conn,spec_list):
        spec_list = '|'.join(spec_list)

        Prelevements_alias = "pr"
        Patients_alias = "p"

        join_column = 'ID_PATIENT'

        lspq_name = "Laboratoire de santé publique du Québec"
        lspq_addr = "20045, chemin Sainte-Marie, Sainte-Anne-de-Bellevue, QC, Canada"

        authors = "Sandrine Moreira, Ioannis Ragoussis, Guillaume Bourque, Jesse Shapiro, Mark Lathrop and Michel Roger on behalf of the CoVSeQ research group (http://covseq.ca/researchgroup)"

        SUBMITTER = "\'SandrineMoreiraLSPQ\' as submitter"
        FN = "\'all_sequences.fasta\' as fn"
        COVV_VIRUS_NAME = "{0}.GENOME_QUEBEC_REQUETE as covv_virus_name".format(Prelevements_alias) #decoration ajouter dans CovGisaidSubmissionV2.py
        COVV_TYPE = "\'betacoronavirus\' as covv_type" 
        COVV_PASSAGE = "\'Original\' as covv_passage"
        COVV_COLLECTION_DATE = "{0}.DATE_PRELEV_HOPITAL as covv_collection_date".format(Prelevements_alias)
        COVV_LOCATION = "\'North-America / Canada / Quebec\' as covv_location"
        COVV_ADD_LOCATION = "\' \' as covv_add_location"
        COVV_HOST = "\'Human\' as covv_host" 
        COVV_ADD_HOST_INFO = "\' \' as covv_add_host_info"
        COVV_GENDER = "{0}.SEXEINFO as covv_gender".format(Patients_alias)

        DTNAISSINFO = "{0}.DTNAISSINFO".format(Patients_alias) # va permettre de determiner l age

        COVV_PATIENT_STATUS = "\'unknown\' as covv_patient_status"
        COVV_SPECIMEN = "\' \' as covv_specimen"
        COVV_OUTBREAK = "\' \' as covv_outbreak"
        COVV_LAST_VACCINATED = "\' \' as covv_last_vaccinated"
        COVV_TREATMENT = "\' \' as covv_treatment"
        
        #covv_seq_technology obtenu dans dans CovGisaidSubmissionV2.py

        COVV_ASSEMBLY_METHOD = "\' \' as covv_assembly_method"
        COVV_COVERAGE = "\' \' as covv_coverage"
        COVV_ORIG_LAB = "{0}.NOM_HOPITAL as covv_orig_lab".format(Prelevements_alias)
        COVV_ORIG_LAB_ADDR = "{0}.ADRESSE_HOPITAL as covv_orig_lab_addr".format(Prelevements_alias)
        COVV_PROVIDER_SAMPLE_ID = "\' \' as covv_provider_sample_id"
        COVV_SUBM_LAB = "'" + lspq_name  + "' as covv_subm_lab"
        COVV_SUBM_LAB_ADDR = "'" + lspq_addr + "' as covv_subm_lab_addr"
        COVV_SUBM_SAMPLE_ID = "{0}.GENOME_QUEBEC_REQUETE as covv_subm_sample_id".format(Prelevements_alias) 
        COVV_AUTHORS = "'" + authors + "' as covv_authors"


        sql = "SELECT {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25} FROM Prelevements {26} inner join Patients {27} on {27}.{28} = {26}.{28} WHERE {26}.GENOME_QUEBEC_REQUETE REGEXP '{29}'".format(SUBMITTER,FN,COVV_VIRUS_NAME,COVV_TYPE,COVV_PASSAGE,COVV_COLLECTION_DATE,COVV_LOCATION,COVV_ADD_LOCATION,COVV_HOST,
                                     COVV_ADD_HOST_INFO,COVV_GENDER,DTNAISSINFO,COVV_PATIENT_STATUS,COVV_SPECIMEN,COVV_OUTBREAK,COVV_LAST_VACCINATED,
                                     COVV_TREATMENT,COVV_ASSEMBLY_METHOD,COVV_COVERAGE,COVV_ORIG_LAB,COVV_ORIG_LAB_ADDR,COVV_PROVIDER_SAMPLE_ID,COVV_SUBM_LAB,
                                     COVV_SUBM_LAB_ADDR,COVV_SUBM_SAMPLE_ID,COVV_AUTHORS,Prelevements_alias,Patients_alias,join_column,spec_list)


        df = pd.read_sql(sql,con=conn)
        #print("DF SHAPE IS ", df.shape)
        '''
        if df.shape[0] == 0:
            print("NO DATA FOR ",spec_list)
        '''
        return df

class Utils():

    @classmethod
    def GetFrenchMonth(cls,english_month):
        month_map = {'January':'Janvier','February':'Février','March':'Mars','April':'Avril','May':'Mai','June':'Juin','Juillet':'July','August':'Août',
                     'September':'Septembre','October':'Octobre','November':'Novembre','December':'Décembre'} 

        return(month_map[english_month])
