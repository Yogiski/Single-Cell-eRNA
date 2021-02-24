import os
import subprocess
import argparse
import pandas as pd

class Controller:

        def __init__(self, *args, out=None):

                self.patient = args[0]
                self.sample_path = args[1]
                self.data_path = args[2]
                self.labeled = args[3]

                if out == None:
                        self.out = os.getcwd()
                else:
                        self.out = out

        def find_atac(self, info):

                samp_info = self.get_tissue_date_id(info, "ATAC")
                atac = self.data_path + "/scATAC-seq_{0}.{1}/{2}/outs/".format(*samp_info)

                return atac

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

                print("Generating file paths")
                samples = pd.read_csv(self.sample_path, index_col="Patient")

                info = samples.loc[self.patient]
                atac = self.find_atac(info)
                rna_bam = self.find_rna_bam(info)
                rds = self.find_labeled_rnaseq()

                print("Calling Label Transfer and Differential Accesibility Script")
                transfer_process = subprocess.Popen(["Rscript", "label-transfer-da.r", atac, rds, self.out])

