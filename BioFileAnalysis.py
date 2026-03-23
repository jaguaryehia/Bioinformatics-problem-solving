from BioHandling import *


class AnalysisFastq:

    def __init__(self, fastq_file):
        self.fastq = read_FASTQ_File(fastq_file)

    def get_Info(self):
        return get_FASTQ_info(self.fastq)
