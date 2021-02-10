import os
import argparse
from controller import Controller

if __name__ == "__main__":

        p = argparse.ArgumentParser(description="""\
                                     Detect eRNAs in scRNA seq data\n
                                     ______________________________\n
                                     \n
                                     """)

        p.add_argument("patient_id",
                       metavar = 'id',
                       type=str,
                       help="Patient ID for file identification")

        p.add_argument("sample_meta",
                       metavar = "s",
                       type=str,
                       help="Path to single cell sample summary file")

        p.add_argument("datastore",
                       metavar = "d",
                       type=str,
                       help="Path to franco lab datastore")

        p.add_argument("scrna_labeled",
                       metavar = "l",
                       type=str,
                       help="Path to labeled scRNA-seq data")

        args = p.parse_args()

        cntrl = Controller(args.patient_id, args.sample_meta, args.datastore, args.scrna_labeled)
        cntrl.run()

        
