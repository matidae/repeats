import sys
import random
import yaml
import numpy
from Bio import SeqIO

class Repeat:
    def __init__(self, name, sequence, num_rep, identity, sd, indels, tsd, frag):
        self.name = name
        self.sequence = sequence
        self.num_rep = num_rep
        self.identity = identity
        self.sd = sd
        self.indels = indels
        self.tsd = tsd
        self.frag = frag

#Load params from YAML file config.yml in same directory
def parse_yaml():
    params = yaml.load(open('config.yml', 'r'))#, Loader=yaml.FullLoader)
    return params

#Load collection of repeats and params
def load_repeats(params):
    repeats_dict = {}
    fasta = SeqIO.index(params['rep_fasta'], "fasta")
    with open(params['rep_list'], 'r') as repeats_file:
        next(repeats_file)
        for line in repeats_file:
            elem = line.rstrip().split()
            name = elem[0]
            sequence = str(fasta[name].seq).upper()
            num_rep = int(elem[1])
            identity = int(elem[2])
            sd = int(elem[3])
            indels = int(elem[4])
            tsd = True if elem[5] == "y" else False
            frag = int(elem[7])
            repeat = Repeat(name,sequence, num_rep, identity, sd, indels, tsd, frag)
            repeats_dict[name] = repeat
    return repeats_dict

###Load other variables###
#Calculate length of sequence of all repeats 
#sum_rep_length = sum([len(rep.sequence) * rep.num_rep for rep in repeats])
#Calculate length of sequence that is going to be randomly generated
#rand_seq_length = seq_length - sum_rep_length

def generate_random_sequence(params):
    #Create DNA alphabet for random sequence
    alphabet = ["T", "G", "C", "A"]
    #Generate random sequence that is going to separate the repeats
    base_sequence = "".join([random.choice(alphabet) for i in xrange(0,params['seq_length'])])
    return base_sequence

#Randomly select positions to insert all the repeats
def assign_coord_repeats(params, repeats_dict):
    total_num_rep = sum([repeats_dict[rep].num_rep for rep in repeats_dict])
    random_repeats_coords = random.sample(xrange(params['seq_length']-1), total_num_rep)
    random_repeats_coords.sort()
    return random_repeats_coords

def shuffle_repeats(repeats_dict):
    #Generate list of all inserts
    total_names_repeats_tmp = [repeats_dict[rep].num_rep*[repeats_dict[rep].name] for rep in repeats_dict]
    #fix: flatten list
    total_names_repeats = [j for i in total_names_repeats_tmp for j in i]
    #Shuffle sequence names.
    random.shuffle(total_names_repeats)
    return total_names_repeats

#Get identity using a normal distribution
def get_identity(mean, sd):
    identity = int(numpy.random.normal(mean, sd, 1))
    while  identity > 100:
        identity = int(numpy.random.normal(mean, sd, 1))
    return identity

#Generate vector of coords for base_changes and indels
def generate_mismatches(sequence, identity, indels):
    alphabet = ["T", "G", "C", "A"]
    seq_len = len(sequence)
    seq = sequence
    num_changes = seq_len - int(round((seq_len * identity/100.0)))
    pos_changes_vec = random.sample(xrange(seq_len), num_changes)
    num_indels = int(round(num_changes * (indels/100.0)))
    indel_changes_vec = random.sample(pos_changes_vec, num_indels)
    base_changes_vec = list(set(pos_changes_vec) - set(indel_changes_vec))
    base_changes_vec.sort()
    indel_changes_vec.sort()
    return base_changes_vec, indel_changes_vec

def add_base_changes(repeat_seq, base_changes_vec):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    for pos in base_changes_vec:
        new_base = random.choice(list(set(alphabet) - set(repeat_seq_list[pos])))
        repeat_seq_list[pos] = new_base
    new_repeat_seq =  "".join(repeat_seq_list)
    return new_repeat_seq

##Add indels
def add_indels(repeat_seq, indels_changes_vec):
    alphabet = ["T", "G", "C", "A"]
    repeat_seq_list = list(repeat_seq)
    for i in xrange(len(indels_changes_vec)):
        if random.choice([0,1]):
            new_base = random.choice(alphabet)
            pos = indels_changes_vec[i]
            repeat_seq_list.insert(pos, new_base)
            for j in xrange(len(indels_changes_vec)):
                indels_changes_vec[j] +=1
        else:
            repeat_seq_list.pop(i)
            for j in xrange(len(indels_changes_vec)):
                indels_changes_vec[j] -=1
    new_repeat_seq =  "".join(repeat_seq_list)
    return new_repeat_seq

def create_TSD(identity, indels):
    alphabet = ["T", "G", "C", "A"]
    tsd_seq_5 = "".join([random.choice(alphabet) for i in xrange(random.choice(xrange(5, 20)))])
    tsd_len = len(tsd_seq_5)
    tsd_base_changes_vec, tsd_indels_changes_vec  = generate_mismatches(tsd_seq_5, identity, indels)
    tsd_seq_mismatches = add_base_changes(tsd_seq_5, tsd_base_changes_vec)
    tsd_seq_3 = add_indels(tsd_seq_mismatches, tsd_indels_changes_vec)
    return tsd_seq_5, tsd_seq_3

def fragment (seq, frag):
    print frag
    frag_size = random.randint(20,95)
    len_seq =len(seq)
    cut = int(len_seq*(frag_size/100.0))
    print cut

##Generate new sequence including the repeats in the random one
def generate_sequence(repeats_dict, rand_rep_pos, rand_seq, total_names_rep):
    seq = ""
    tsd_seq_5= ""
    tsd_seq_3= ""
    pre_n = 0
    n=0
    new_repeats_coord = []
    for n,m in zip(rand_rep_pos, total_names_rep):
        #Get sequence of repeat
        repeat_seq = repeats_dict[m].sequence
        #Get identity from a normal distribution
        identity = get_identity(repeats_dict[m].identity, repeats_dict[m].sd)
        #Get base_changes and indels vectors and identity
        identity_fix = identity + (100 - identity) * 0.5
        print identity, m
        base_changes_vec, indels_changes_vec = generate_mismatches(repeats_dict[m].sequence, identity_fix, repeats_dict[m].indels)
        #Add mismatches to original repeat sequence
        repeat_seq_mismatches = add_base_changes(repeat_seq, base_changes_vec)
        #Add indels to original repeat sequence
        new_repeat_seq = add_indels(repeat_seq_mismatches, indels_changes_vec)
        new_repeat_seq = new_repeat_seq.lower()
        #Check if TE creates TSDs
        if repeats_dict[m].tsd:
            #Generate TSD
            tsd_seq_5, tsd_seq_3 = create_TSD(identity_fix, repeats_dict[m].indels)
        #Append new repeat sequence to base sequence
        seq += rand_seq[pre_n:n] + tsd_seq_5 + new_repeat_seq + tsd_seq_3
        #Get new repeat sequence end coordinate
        repeat_end = len(seq) - len(tsd_seq_3)
        #Get new repeat sequence start coordinate
        repeat_start = repeat_end - len(new_repeat_seq) + 1
        #---------------
        new_repeat_seq = fragment(new_repeat_seq, repeats_dict[m].frag)

        #Append to vector new data about new repeat useful for a GFF
        new_repeats_coord.append([str(repeat_start), str(repeat_end), new_repeat_seq, identity])
        #Sets new end coordinate as start for next roung
        pre_n = n
    #At the last step add the remaining base sequence
    seq += rand_seq[n:]
    #Return sequences and repeat data
    return seq, new_repeats_coord

#Print final sequence to stdout
def print_data(prefix, seq, new_repeats_coord, total_names_rep):
    fasta_out = open(prefix + "_out_sequence.fasta", "w")
    fasta_out.write( ">sequence\n" )
    for n in xrange(0,len(seq),100):
        fasta_out.write(str(seq[n:n+100]) + "\n")
    fasta_out.close()
    #Print start positions of repeats to stderr
    fasta_rep = open(prefix + "_out_repeats.fasta", "w")
    gff_rep = open(prefix + "_out_repeats.gff", "w")
    c = 1
    for n,m in zip(new_repeats_coord, total_names_rep):
        #rep_name = ">" + m + "_" +  n[0] + "_" + n[1] + "\n"
        repeat_identity = str(n[3])
        repeat_name = ">" + m + "_p" +  str(c) +  "_" + repeat_identity + "\n"
        repeat_sequence = str(n[2])+ "\n"

        fasta_rep.write(repeat_name)
        fasta_rep.write(repeat_sequence)

        gff_rep.write("\t".join(["sequence", "script", "repeat_region", n[0], n[1], ".", ".", ".", "ID=" + m + "_p" + str(c) + ";identity=" + repeat_identity + "\n"]))
        c+=1
    fasta_rep.close()
    gff_rep.close()

def main():
    #Load parameters
    params = parse_yaml()
    #Load repeat sequences
    repeats_dict = load_repeats(params)
    #Generate random sequence
    base_sequence = generate_random_sequence(params)
    #Assign random positions to repeats
    repeats_coord = assign_coord_repeats(params, repeats_dict)
    #Shuffle repeats do they aren't  in order 
    shuffled_repeats = shuffle_repeats(repeats_dict)
    #Insert repeats in sequence (after applying identity changes), report new sequence and coordinates
    sequence, new_repeats_coord = generate_sequence(repeats_dict, repeats_coord, base_sequence, shuffled_repeats)
    #Output fasta file with new sequence, repeats, and GFF.
    print_data(params['prefix'], sequence, new_repeats_coord, shuffled_repeats)

if __name__ == "__main__":
    main()
