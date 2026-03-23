import argparse
import datetime
from pathlib import Path
from BioIO import *


class Args:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description="This a bioinformatic command line tool help into the biological data analysis")

        self.parser.add_argument("-file", "--path", help="return file that has known type functions")
        self.parser.add_argument("-D", "--dna", help="return dna functions", type=str)
        self.parser.add_argument("-R", "--rna", help="return rna functions", type=str)
        self.parser.add_argument("-p", "--protein", help="return protein functions", type=str)
        # self.parser.add_argument("-tcript","--transcript",help="return transcript of dna or rna sequence", type=str)
        # self.parser.add_argument("-tlate","--translation",help="return translate of dna or rna sequence", type=str)
        # self.parser.add_argument("-comp","--complement",help="return complement", type=str)
        # self.parser.add_argument("-rcomp","--reversecomplement",help="return Reverse Complement", type=str)
        # self.parser.add_argument("-GC","--GCcontent",help="return GC Content", type=str)
        # self.parser.add_argument("-btcript","--backtranscript",help="return back Transcript", type=str)
        self.args = self.parser.parse_args()
        self.bio = None

    def getArgs(self):
        return self.args

    def printFunction(self, data):
        if data == 'dna':
            print(
                '1. Translation\n2. Transcription\n3. Complement\n4. Reverse Complement\n5. GC Content\n6. Back Transcript')

            self.parser.add_argument("-tcript", "--transcript", help="return transcript of dna or rna sequence",
                                     type=str)
            self.parser.add_argument("-tlate", "--translation", help="return translate of dna or rna sequence",
                                     type=str)
            self.parser.add_argument("-comp","--complement",help="return complement", type=str)
            self.parser.add_argument("-rcomp","--reversecomplement",help="return Reverse Complement", type=str)
            self.parser.add_argument("-GC","--GCcontent",help="return GC Content", type=str)
            self.parser.add_argument("-btcript","--backtranscript",help="return back Transcript", type=str)
        elif data == 'rna':
            print('rna data')
            '''
            translation
            rcomplement
            complement
            back transcript
            
            '''
        elif data == 'protein':
            '''complement
                rcomplement
                 
            '''
            print('pro data')
        else:
            print(f'What do you wanna do to your file: ')
        # print('7. transcription')
        # print('8. transcription')
        # print('9. transcription')
        # print('10. transcription')

    def DataType(self):
        if self.args.dna:
            print(f'What function do you wanna do on this seq: {self.args.dna}')
            self.bio = BioIO(self.args.dna)
            self.printFunction('dna')

        elif self.args.rna:
            print(f'What function do you wanna do on this seq: {self.args.rna}')
            self.printFunction('rna')
        elif self.args.protein:
            print(f'What function do you wanna do on this seq: {self.args.protein}')
            self.printFunction('protein')
        elif self.args.file:
            print(f'What function do you wanna do on this file: {self.args.file}')
            self.printFunction('file')
        else:
            print(f'Wrong arguments')
