from Bio.Blast import NCBIWWW
import pandas as pd
import numpy as np
import textwrap
from Bio.Blast import NCBIXML
import socket
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
from Bio import Entrez
from Bio.Seq import Seq

class entrez:
    def __init__(self,target,email_for_errors):
        self.target = target
        self.records = None
        self.email = email_for_errors
        self.seq = None
    def entrez_get_seq(self):
        Entrez.email = self.email
        handle = Entrez.efetch(db="nucleotide", id=self.target, retmode="xml",rettype="gb")
        records = Entrez.parse(handle)
        for record in records:
            self.seq = str(Seq(record['GBSeq_sequence']))
        handle.close()




class blast:

    def __init__(self, fname, blast = False): # Input file name as acdefgh for abcdefgh.xml or .fasta
        self.fname_in = fname + ".xlsx"
        self.fname_fasta = fname + ".fasta"
        self.fname_results = fname + "results.xml"
        self.data = None
        self.result_handle = None
        self.test_seq = []
        self.output = None
        self.load_data()
        self.split_data()
        if blast == True:
            self.write_fasta()
            self.blast_seq()
            self.save_results()

    def load_data(self):
     # Load Data
        data = pd.read_excel(self.fname_in, engine='openpyxl')
        self.data = data.to_numpy()
        self.output = np.copy(self.data)
        self.raw_seq = self.data[2:, 1]
        self.probe = self.data[2:, 0]

    def split_data(self):

        for seq in self.raw_seq:
            self.test_seq.append((seq.split())[-1])

    def write_fasta(self):

        with open(self.fname_fasta, 'w') as f:
            for i, seq in enumerate(self.test_seq):
                f.write(f">{self.probe[i]}\n")
                f.write(f"{seq}\n\n")

    def blast_seq(self):

        sequence_data = open(self.fname_fasta).read()
        self.result_handle = NCBIWWW.qblast(
            "blastn", "GPIPE/9606/current/ref_top_level GPIPE/9606/current/rna", sequence_data, megablast=True)

    def save_results(self):
        with open(self.fname_results, 'w') as save_file:
            blast_results = self.result_handle.read()
            save_file.write(blast_results)


class blast_result:

    def __init__(self, fname_results, matches = None):
        self.fname_results = fname_results
        self.relevant_results = []
        self.matches = matches
        self.check_relevant_results(self.fname_results,self.matches)
        self.print_all_results()

    def print_all_results(self):
        E_VALUE_THRESH = 100
        for record in NCBIXML.parse(open(self.fname_results)):
            if record.alignments:
                print("\n")
                print("query: %s" % record.query[:100])
                for align in record.alignments:
                    for hsp in align.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            print("match: %s " % align.title[:100])

    def check_relevant_results(self, fname_results, matches):

        for i, record in enumerate(NCBIXML.parse(open(fname_results))):
            self.relevant_results.append(True)
            for align in record.alignments:
                if align.accession not in matches:
                    self.relevant_results[i] = False


def calculate_GC_and_mt(seqeunces):
    GC_percent = []
    mt_result = []
    for seq in seqeunces:
        GC_percent.append(GC(seq.replace(" ","")))
        mt_result.append(mt.Tm_GC( seq.replace(" "," "), Na=50))
    return GC_percent,mt_result

# class test_sequences:
#     def __init__(seq,nuc):
#         self.full_seq = seq
#         self.nuc = nuc
#     def split_sequences(self, nuc = self.nuc):
#         for i in range(len(seq)-nuc+1):



    # blast_records[0].alignments[1].accession


if __name__ == '__main__':

    # blast_1 = blast(sys.argv[1])
    # blast_1_results = blast_result(sys.argv[1]+"results.xml",['NR_122046','NC_000022','NC_000002'])
    # blast_results = blast_1_results.relevant_results
    # GC_percent,mt_result = calculate_GC_and_mt(blast_1.raw_seq.tolist())
    # print(mt_result)

    ent = entrez('NR_122046.1','pranay99jain@gmail.com')
    ent.entrez_get_seq()
    print(ent.seq)

#   print(data['Complete Probes (5\'-3\')'])
