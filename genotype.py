"""
Written by Alexander James Kern on May 8th, 2019
Multiplex Genotyping coding challenge by Myriad Genetics

One kind of test Counsyl performs is known as genotyping; genotyping can tell you if a given patient sample has the
normal or the mutant (disease-causing) sequence at one point in the genome. Some kinds of genotyping chemistries have an
interesting property: if you mix up many samples and genotype the mixture the test will read "mutant" if ANY sample in
the mix was mutant, and "NORMAL" only if ALL samples were normal.

For most diseases, the mutant allele (version of the gene) is extremely rare - most samples will test as normal. This
suggests an interesting efficiency trade-off. If you pool combinations of samples in the right way (a strategy known as
"multiplexing"), you can perform fewer tests than the number of samples and use software to find out which individual
sample(s) had a mutation and which were normal. In this problem, you'll write the code that does that computation.

Problem
You will be a given a series of test cases, where each test case describes a series of mixture experiments. Each
mixture experiment will contain 1 or more samples and a reported genotype (NORM or MUT) for that mixture. Following the
multiplex genotyping rule described above, compute the genotype for each sample. Beware! Sometimes the chemistry will
go awry in a batch of samples, and you'll get a result that doesn't make sense. If there is no sample-to-genotype
mapping that is consistent with the multiplexing rule, report an error. It is also possible that the multiplexing was
done poorly, and it's impossible to compute a unique sample-to-genotype map that is consistent with the rules; in this
case you must also report an error.

Input
The input file has several separate test cases (sets of mixtures and their genotype), each terminated by a single
empty line. Treat each test case separately. For example, sample IDs might repeat across separate test cases, but they
are unrelated. Your program can either read the input file directly or accept the file piped in on standard input.

Within a test case, each line (representing a mixture) will be a comma-separated set of fields in the following format:

genotype,sample_id_1,sample_id_2,...,sample_id_n
genotype represents the result of genotyping the mixture and will either be MUT for mutant or NORM for normal. The
sample IDs that follow represent the samples in that particular mixture. They will be non-negative integers. There
will be at least one sample ID per mixture.

Output
There are three possible cases:

1. If a unique mapping of samples-to-genotype exists, output the number of mutant and normal samples followed by the
sample-to-genotype mapping (in sorted order):

MUT COUNT: 2
NORM COUNT: 3
23,MUT
204,NORM
208,NORM
340,MUT
344,NORM
2. If no consistent sample-to-genotype map is possible, output a single line containing only the string INCONSISTENT

3. If more than one map would be consistent, output a single line containing only the string NONUNIQUE.

Terminate the results of each test case with a single empty line.

Sample Input:

NORM,0,1
NORM,1,2
NORM,0,2

NORM,100,110
MUT,110,12

NORM,0,1
MUT,1,2
NORM,1,3
NORM,2,3

MUT,0,1
MUT,1,2

Sample Output:

MUT COUNT: 0
NORM COUNT: 3
0,NORM
1,NORM
2,NORM

MUT COUNT: 1
NORM COUNT: 2
12,MUT
100,NORM
110,NORM

INCONSISTENT

NONUNIQUE
"""
import doctest


def breakFileIntoInputs(file_name):
    """
    Gets file, breaks up each input individually
    >>> inputs = breakFileIntoInputs("input.txt")
    >>> print(inputs[0])
    NORM,80,74
    NORM,32,74,19
    NORM,32,19,75
    NORM,12,28,54
    NORM,74,54,54,12
    NORM,55,28,54,32
    NORM,75,28,96,74
    NORM,55,28,32,80,32
    NORM,12,32,75
    NORM,12,28
    """
    file = open(file_name, "r")
    string_of_inputs = file.read()
    norm_and_mut_inputs = string_of_inputs.split('\n\n')

    # last line appears to be empty, only take 7 of the elements
    return [i for i in norm_and_mut_inputs if len(i) != 0]


def  mapNormAndMut(norm_and_mute_lines):
    """
    note that any norm sequecces found on an input are assumed
    :param norm_and_mute_lines: This is a string of the form:
    "NORM,a1,a2,...,an
     MUT,b1,b2,...,bk
     :
     '
     NORM,z1,z2,..,zj"
    Each line is of an inconsistent size, and is either NORM or MUT leading it.
    For NORM, we assume that there is no mutated sequince id's thus we concatenate all
    of them in the first run. this function takes O(n) time (each iteration of the for
    loop takes a single line, then breakes them up by commas which is O(n) time, than looks
    at each sample which is O(n) time). The nested for loops operate specifically once on
    each line, no more
    :return:
    .
    .
    .
    >>> inputs = breakFileIntoInputs("input.txt")
    >>> x, y = mapNormAndMut(inputs[3])
    >>> print(y)
    [[174], [], [174]]
    >>> print(x)
    {8: 'NORM', 23: 'NORM', 113: 'NORM', 194: 'NORM', 120: 'NORM', 187: 'NORM', 5: 'NORM', 20: 'NORM', 82: 'NORM', 21: 'NORM', 66: 'NORM', 118: 'NORM', 136: 'NORM', 182: 'NORM', 87: 'NORM', 122: 'NORM', 154: 'NORM', 19: 'NORM', 174: 'MUT'}
    """
    norm_map = {}
    norms_and_muts = norm_and_mute_lines.split('\n')
    mutate_lists = []

    # maps all normal functions to norm,
    for line in norms_and_muts:
        title_and_samples = line.split(',')
        if title_and_samples[0] == 'NORM':
            for sample in title_and_samples[1:]:
                # we call samples we KNOW are normal genes
                norm_map[int(sample)] = 'NORM'

    # if it is in a MUT case, make sure it is not in norm
    for line in norms_and_muts:
        title_and_samples = line.split(',')
        if title_and_samples[0] == 'MUT':
            mutate_sequence = []
            for sample in title_and_samples[1:]:
                sequence_id = int(sample)
                if sequence_id in norm_map and norm_map[sequence_id] == 'NORM':
                    continue
                else:

                    # For later use, if this case works we have the list ready to run
                    norm_map[sequence_id] = 'MUT'
                    mutate_sequence.append(sequence_id)
            mutate_lists.append(mutate_sequence)

    return norm_map, mutate_lists


def validMutationList(mutated_lists):
    """
    :param mutated_lists: A list of lists, containg assumed valid mutation
    :return: True if a valid mutation
        The nature of this is each set is a bipartate graph from id to the set it contains, and
    from the set to the ids it contains.
        If a set of ids has n unique ids, in n or more unique sets, then it is the case
    that each id can be uniquely identified as mutated, else we can not make the assumption
    all the ids are mutated, and it must return false.
    Since this uses depth first search, let V be the number of sets plus the number of unique ids.
    The edges are at most (V/2)^2 for a bitpartate graph, but that is usually unlikely, thus the runtime
    is O(V + E)
    .
    .
    .
    >>> validMutationList([[1, 2], [1, 3]])
    False
    >>> validMutationList([[1, 2], [1, 2]])
    True
    >>> validMutationList([[0,1], [1, 2], [2, 1], [3, 4, 5], [3, 4, 5]])
    False
    >>> validMutationList([[0,1], [1, 2], [2, 1], [3, 4, 5], [3, 4, 5], [3, 4, 5]])
    True
    """
    ids_to_sets = {}
    visited_nodes = {}

    # This for loop creates the graph to perform the search. it chunks out the lists in the outer
    # loop, than works on thos chunks in the inner loop, thus it is O(n) runtime
    for i in range(len(mutated_lists)):
        temp_set = mutated_lists[i]
        ids_to_sets[i*-1 - 1] = set()
        for id in temp_set:
            if id not in ids_to_sets:
                ids_to_sets[id] = {i*-1 - 1}
                ids_to_sets[i*-1 - 1].add(id)
                visited_nodes[i*-1 - 1] = False
                visited_nodes[id] = False
            else:
                ids_to_sets[id].add(i* -1 - 1)
                ids_to_sets[i*-1 - 1].add(id)
                visited_nodes[i*-1 - 1] = False
                visited_nodes[id] = False
    for unique_id in visited_nodes:
        if visited_nodes[unique_id]:
            continue
        if not depthFirstSearch(unique_id, ids_to_sets, visited_nodes):
            return False
    return True


def depthFirstSearch(current_id, ids_to_sets, visited_nodes):
    """
    Written from scratch, as this was a bit of a unique situation
    Basically crawls throu ids to sets, if the current id is negative, it is a
    set to a id, if postiive the other. Basic recursive definition.
    O(V + E) time
    :param current_id: current id we are to iterate on
    :param ids_to_sets:
    :param visited_nodes:
    :return: True if # of sets >= # unique ids in particular bipartate graph, else false
    .
    .
    .
    >>> # case where same ids and sets
    >>> ids_to_sets = {1:{-1, -2}, -1: {1, 2}, 2: {-1, -3}, -2: {1 ,3}, 3: {-2, -3}, -3: {3, 2}}
    >>> visited_nodes = {1: False, -1: False, 2: False, -2: False, 3: False, -3: False }
    >>> depthFirstSearch(1, ids_to_sets, visited_nodes)
    True
    >>> # case where more ids than sets
    >>> ids_to_sets = {1:{-1, -2}, -1: {1, 2}, 2: {-1}, -2: {1 ,3}, 3: {-2}, }
    >>> visited_nodes = {1: False, -1: False, 2: False, -2: False, 3: False, -3: False }
    >>> depthFirstSearch(1, ids_to_sets, visited_nodes)
    False
    >>> # case where less ids than sets
    >>> ids_to_sets = {1:{-1, -2, -4}, -1: {1, 2}, 2: {-1, -3}, -2: {1 ,3}, 3: {-2, -3}, -3: {3, 2}, -4:{1, 2, 3}}
    >>> visited_nodes = {1: False, -1: False, 2: False, -2: False, 3: False, -3: False , -4: False}
    >>> depthFirstSearch(1, ids_to_sets, visited_nodes)
    True
    """

    # Keeps track of number of unique sets visited and number of uniqe ids visited.
    def depthHelper(current_id, ids_to_sets, visited_nodes, count):
        if visited_nodes[current_id]:
            return
        if current_id < 0:
            count[0] += 1
        else:
            count[1] += 1
        visited_nodes[current_id] = True
        for unique_id in ids_to_sets[current_id]:
            depthHelper(unique_id, ids_to_sets, visited_nodes, count)
    count = [0, 0]
    depthHelper(current_id, ids_to_sets,  visited_nodes, count)

    # if there are more sets (or equal) than unique items in the bipartate graph,
    # than each id must be a mutated id
    return count[0] >= count[1]


def getOutput(file_name):
    """
    Simply obtains output
    :param file_name:
    :return: answers of file
    """
    inputs = breakFileIntoInputs(file_name)
    output = []
    for input_ids in inputs:
        x, y = mapNormAndMut(input_ids)
        if len(y) == 0:
            output.append(x)
        elif [] in y:
            output.append('INCONSISTENT')
        elif validMutationList(y):
            output.append(x)
        else:
            output.append('NONUNIQUE')

    return output


def writeOutput(output, file_name):
    """
    :param output: results of the Multiplex Genotyping
    :param file_name: name of file
    :return: None
    """
    # This is just some standard output file tequnqes applied to the solution. Since the keys need to
    # be sorted it runs in O(nlogn) time.
    f = open(file_name, "w")
    for result in output:
        if result == 'INCONSISTENT':
            f.write('INCONSISTENT\n')
        elif result == 'NONUNIQUE':
            f.write('NONUNIQUE\n')
        else:
            sorted_keys = list(result.keys())
            sorted_keys.sort()
            num_norms = 0
            num_mutes = 0
            for key in sorted_keys:
                if result[key] == 'NORM':
                    num_norms += 1
                else:
                    num_mutes += 1
            f.write("MUT Count: {}\n".format(num_mutes))
            f.write("NORM Count: {}\n".format(num_norms))
            for key in sorted_keys:
                f.write(str(key))
                f.write(',')
                f.write(result[key])
                f.write('\n')
        f.write("\n")
    f.close()

    # Remove last two empty lines
    read_file = open(file_name)
    lines = read_file.readlines()
    read_file.close()
    w = open(file_name, 'w')
    w.writelines([item for item in lines[:len(lines) - 2]])
    lastline = lines[len(lines) - 2]
    w.write(lastline[:len(lastline) - 1])
    w.close()
    return None


if __name__ == "__main__":
    out = getOutput("input.txt")
    writeOutput(out, "output.txt")
    doctest.testmod()
