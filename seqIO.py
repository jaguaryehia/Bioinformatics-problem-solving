import re

'''
    function for reading HUMMER output file return list of id desired using regex
        Input:
            :param input : filepath 
        Output:
            :return list() : list of ids to be used again
'''


def readHUMMEROutput(filepath):
    pattern = re.compile(r"\w+_\d*.\d")
    id_list = []
    with open(filepath, 'r') as f:
        for line in f:
            id_found = pattern.findall(line)
            if id_found:
                for id in id_found:
                    id_list.append(id)
    f.close()
    return id_list


'''
    function for writing the output file after processes 
        Input:
            :param dict(): dictionary storing the ids and seqs
            :param outputfile: output file name
            :param formatter: formatter for outputfile
        Output:
            :return file: output file that have the seqs and ids as fasta file 
'''


def writeSeqs(seq_dict, outputfile, formatter):
    with open(outputfile, 'w') as output:
        for id_seq in seq_dict:
            seq = seq_dict[id_seq]
            seq = [seq[i:i + formatter] for i in range(0, len(seq), formatter)]
            output.write(f">{id_seq}\n")
            for line in seq:
                output.write(f"{line}\n")
    output.close()


def getGeneByIndex(seq_dict, id, start, end):
    return seq_dict[id][start - 1:end]
