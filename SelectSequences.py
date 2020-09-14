import datetime
from Bio import SeqIO
import glob
import os
import re
import shutil

_debug = False

published_seqid = []
published_path = "/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/FinalPublished/"

for published_fasta in glob.glob(published_path + "*/all_sequences.fasta"):
    for rec in SeqIO.parse(published_fasta,'fasta'):
        rec_id = re.search(r'^hCoV-19/Canada/Qc-(L\S+)/\d{4}$',rec.id).group(1)
        published_seqid.append(rec_id)


if _debug:
    source_seq_path = "/data/Applications/GitScript/GisaidSubmission/TEST2/SOURCE/"
    nextstrain_seq = "/data/Applications/GitScript/GisaidSubmission/TEST2/NEXTSTRAIN/subsampled_alignment_quebec.fasta"
    dest_path = "/data/Applications/GitScript/GisaidSubmission/TEST2/DEST/"
else:
    source_seq_path = "/data/Runs/SARS-CoV-2/GenomeCenterSeq/2020_07_06-full/"
    nextstrain_seq = "/data/PROJETS/COVID_19/Freeze1_wExtra_large_no_rta_sampling_corrected/results/subsampled_alignment_quebec.fasta"
    dest_path = "/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease/FinalUnpublished/"

compile_obj = re.compile(r'^Canada/Qc-(\S+)/\d{4}$')

print("Start")

nb_transfered = 0

for myrec in SeqIO.parse(nextstrain_seq,"fasta"):
    search_obj = compile_obj.search(myrec.id)
    if search_obj:
        lspq_id = search_obj.group(1)
        
        if lspq_id not in published_seqid:
            for fasta in os.listdir(source_seq_path):
                if re.search(lspq_id + r'\.consensus\S+(pass|flag).fasta',fasta):
                    nb_transfered += 1
                    shutil.copyfile(os.path.join(source_seq_path,fasta),os.path.join(dest_path,fasta))

print("End transfer ", nb_transfered , " from  ", source_seq_path, " to ", dest_path)
