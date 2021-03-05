import os
import sys
import subprocess
import pandas as pd

class Controller:

    """
    Controller class manages the detect eRNAs pipeline.
    It generates file paths based on user input and sample metadata and verifies their existance.
    It then calls the the various analysis scripts in the pipeline and saves their outputs.
    """

    def __init__(self, *args, out=None):

        self.patient = args[0]
        self.sample_file = args[1]
        self.data_path = args[2]

        if out == None:
                self.out = os.getcwd() + "/detect-ernas-out/"
        else:
                   self.out = out

    def find_atac(self, info):

        samp_info = self.get_tissue_date_id(info, "ATAC")
        atac = self.data_path + "/scATAC-seq_{0}.{1}/{2}/outs/".format(*samp_info)
        files = os.listdir(atac)
        assert "fragments.tsv.gz" in files
        assert "filtered_peak_bc_matrix.h5" in files
        assert "singlecell.csv" in files

        return atac

    def find_rna_bam(self, info):

        samp_info = self.get_tissue_date_id(info, "RNA")
        source = self.data_path + "/scRNA-seq_{0}.{1}/{2}/outs/".format(*samp_info)
        for f in os.listdir(source):
            if ("possorted_" in f) and (".bam" in f) and (".bai" not in f):
                bam = source + f
                return bam

        print("Could not find necessary scRNA-seq possorted bam file. \nExitting..")
        sys.exit()


    def get_tissue_date_id(self, info, style):

        tissue = info["Tissue"]
        date = info["date_{}".format(style)]
        samp_id = info["sampleid_{}".format(style)]

        return tissue, date, samp_id

    def find_labeled_rnaseq(self):

        rds_dir = self.data_path + "/scRNAlabeled-rds-files/"
        files = os.listdir(rds_dir)
        for f in files:
                if self.patient in f:
                        return rds_dir + f

        print("Could not find labeled scRNA rds file for patient {}\nExitting...".format(self.patient))
        sys.exit()

    def make_dir(self, dirstring):

        try:
            os.mkdir(dirstring)
        except FileExistsError:

            if len(os.listdir(dirstring)) != 0:
                usr = input("""Directory {} exists and is not empty.
                Continuing this script may overwrite previously saved files.
                Do you want to proceed? y/n> """.format(dirstring))

                if (usr == "y") or (usr == "Y") or (usr == 'yes') or (usr == "Yes"):
                    pass
                else:
                    sys.exit()
            else:
                pass

    def verify_working_dir(self):

        required = ["label-da-bed.r",
                    "hg38.rds",
                    "GRCh38-ccREs.dELS.bed",
                    "GRCh38-ccREs.pELS.bed",
                    "cellranger-dna-1.1.0",
                    "refdata-GRCh38-1.0.0"]

        files = os.listdir()
        count = 0

        for r in required:
            try:
                assert r in files
                count += 1
            except AssertionError:
                break

        if "Single-Cell-eRNA" in files:
            os.chdir("Single-Cell-eRNA")
            self.verify_working_dir()
        elif count == len(required):
           self.cwd = os.getcwd()
           return
        else:
           print("Could not find files necessary for analysis. Exiting...")
           sys.exit()


    def setup(self):

        print("Gathering Sample Information")
        samples = pd.read_csv(self.data_path + self.sample_file, index_col="Patient")
        info = samples.loc[self.patient]

        print("Generating file paths")
        self.atac = self.find_atac(info)
        self.rna_bam = self.find_rna_bam(info)
        self.rds = self.find_labeled_rnaseq()

        print("Making output directories")
        try:
            os.mkdir(self.out)
        except FileExistsError:
            pass

        self.make_dir(self.out + "/" + self.patient)

        self.bed_dir = "{0}/{1}/da-enhancer-beds".format(self.out, self.patient)
        self.make_dir(self.bed_dir)

        self.slice_dir = "{0}/{1}/sliced-rna-bams".format(self.out, self.patient)
        self.make_dir(self.slice_dir)

        self.ids_dir = "{0}/{1}/rna-celltype-ids".format(self.out, self.patient)
        self.make_dir(self.ids_dir)

        self.erna_dir = "{0}/{1}/erna-coverage".format(self.out, self.patient)
        self.make_dir(self.erna_dir)

    def run_analysis(self):
        print("Starting Label Transfer, DA Analysis, and BED file generation\nThis may take a while...")
        subprocess.run(["Rscript", "label-da-bed.r", self.atac, self.rds, self.out, self.patient])

    def run_bamslice(self):

        env = os.environ.copy()
        env["PATH"] = "{0}/cellranger-dna-1.1.0:{1}".format(self.cwd, os.environ["PATH"])
        command = "cellranger-dna bamslice --id {0} --csv {1} --bam {2}".format("test_subset", self.ids_dir + "/tumor_subsets.csv", self.rna_bam)
        subprocess.Popen(command, shell=True, env=env).wait()

    def run_coverage(self):
        env = os.environ.copy()
        clusters = [i for i in os.listdir(self.cwd + "/test_subset/outs/subsets/") if ".bai" not in i]
        for c in clusters:
            name = c.split(".")[0]
            print("Looking for ernas in {}".format(name))
            a = self.bed_dir + "/" + name + "-da-enhancers.bed"
            b = self.cwd + "/test_subset/outs/subsets/" + c
            out = self.erna_dir + "/" + name + ".txt"
            command = "bedtools coverage -a {0} -b {1} > {2}".format(a, b, out)
            subprocess.Popen(command, shell=True, env=env).wait()

