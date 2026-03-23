from BioHandling import *


class BioIO:
    def __init__(self, data):
        self.Data = data

    def getType(self):
        return seq_Type(self.Data)

    # def dna(self):
