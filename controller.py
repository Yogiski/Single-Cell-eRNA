import os
import argparse
import pandas as pd

class Controller:

        def __init__(self, *args):# sample_path=None, patient, data_path, labeled_path, outpath=None, log=True):

                self.patient = args[0]
                self.sample_path = args[1]
                self.data_path = args[2]
                self.labeled = args[3]

        def find_fragments(self, info):

                samp_info = self.get_tissue_date_id(info, "ATAC")
                fragments = self.data_path + "/scATAC-seq_{0}.{1}/{2}/outs/fragments.tsv.gz".format(*samp_info)

                return fragments

        def find_rna_bam(self, info):

                samp_info = self.get_tissue_date_id(info, "RNA")
                bam = self.data_path + "/scRNA-seq_{0}.{1}/{2}/outs/possorted_bam.bam".format(*samp_info)

                return bam

        def get_tissue_date_id(self, info, style):

                tissue = info["Tissue"]
                date = info["date_{}".format(style)]
                samp_id = info["sampleid_{}".format(style)]

                return tissue, date, samp_id

        def find_labeled_rnaseq(self):

                files = os.listdir(self.labeled)
                for f in files:
                        if self.patient in f:
                                return self.labeled + "/" + f

        def run(self):

                samples = pd.read_csv(self.sample_path, index_col="Patient")

                info = samples.loc[self.patient]
                fragments = self.find_fragments(info)
                rna_bam = self.find_rna_bam(info)
                rds = self.find_labeled_rnaseq()

                print(os.listdir(self.data_path))
                print(fragments)
                print(rds)
                print(rna_bam)

