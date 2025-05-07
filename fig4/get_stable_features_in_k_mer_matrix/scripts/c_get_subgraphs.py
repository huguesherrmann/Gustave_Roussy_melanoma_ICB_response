import networkx as nx
import itertools
import argparse


# ......................................................
#
#   FUNCTIONS ----
#
# ......................................................
def subgraph_counts(lengths_dict, threshold = 1000):
    """
    Given a dictionary of node names and values, this function builds a graph where an edge between two nodes exists
    if the absolute difference of their values is less than a given threshold. It then finds all connected subgraphs
    in the graph, sorts them in ascending order of their minimum node value, and assigns a unique identifier to each
    subgraph. Finally, it returns a dictionary mapping node names to their respective subgraph identifiers.
    
    Parameters:
    lengths_dict (dict): A dictionary mapping node names to values.
    threshold (int, optional): The maximum absolute difference between two node values for an edge to exist. Default is 1000.

    Returns:
    res_dict (dict): A dictionary mapping node names to their respective subgraph identifiers. Nodes in the same subgraph share the same identifier.
    """

    G = nx.Graph()
    
    # Add nodes to the graph with custom names and values
    for name, value in lengths_dict.items():
        # if value != 0:  # Only add the node if its value is not 0
        G.add_node(name, value=value)

    # Add edges between nodes based on the absolute difference of their values
    for node1, node2 in itertools.combinations(G.nodes, 2):
        value1 = G.nodes[node1]['value']
        value2 = G.nodes[node2]['value']
        
        if abs(value1 - value2) < threshold:
            G.add_edge(node1, node2)

    # Create a list of subgraphs, where each subgraph is a connected component of the graph G then sort them by the minimum value of their nodes
    S = sorted([G.subgraph(c).copy() for c in nx.connected_components(G)], key=lambda x: min(data['value'] for node, data in x.nodes(data=True)))
    subgraph_id = 0
    res_dict = dict()
    for sub in S:

        for name,data in sub.nodes(data=True):

            res_dict.update({name:subgraph_id})

        subgraph_id += 1

    
    return res_dict


def parse_file_to_dict(file_path, chromosome):
    result_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into chromosome, key, and value
            chr_, key, value = line.strip().split()
            # Check if the chromosome matches the specified chromosome
            if chr_ == chromosome:
                # Add key-value pair to the dictionary
                result_dict[key] = int(value)
    return result_dict


def write_dict_to_file(dictionary, file_path, prefix = ""):
    with open(file_path, 'a') as file:
        for key, value in dictionary.items():
            file.write(f"{prefix} {key} {value}\n")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser = argparse.ArgumentParser(description = "Draw edges between nodes (contigs) if their distance is less than a threshold.")
parser.add_argument('-i', '--input', type = str, help = 'Input path. A file with tag and genomic position.')
parser.add_argument('-o', '--output', type = str, help = 'Output filename.')
parser.add_argument('-d', '--distance', type = int, default = 1000, help = 'Distance in base pair to connect 2 nodes.')
args = parser.parse_args()
#input_file = "/mnt/beegfs/scratch/h_herrmann/DEV/get_subgraphs/to_connect.txt"
#output_file = '/mnt/beegfs/scratch/h_herrmann/DEV/get_subgraphs/connected.txt'


# ......................................................
#
#   GET SUBGRAPHS ----
#
# ......................................................
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
               "chr21", "chr22", "chrX", "chrY"]
for chromosome in chromosomes:
    result_dict = parse_file_to_dict(args.input, chromosome)
    graph_results = subgraph_counts(result_dict, threshold = args.distance)
    
    write_dict_to_file(graph_results, args.output, prefix = chromosome)
