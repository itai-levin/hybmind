import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
import json
from queue import PriorityQueue
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles
from graph_tool.all import *
import networkx as nx
import glob
import re
import matplotlib as mpl


def show_reaction_ls (reaction_ls, metadata = None):
    for i, reaction in enumerate(reaction_ls):
        try:
            if metadata:
                print (metadata[i])
            plt.figure(figsize=(10,10))
            plt.imshow(Chem.Draw.ReactionToImage(Chem.ReactionFromSmarts(reaction, useSmiles=True), subImgSize=(500,500)))
            plt.xticks([])
            plt.yticks([])
            plt.show()
        except:
            print (reaction)

def get_all_children(tree, field = 'smiles', reactions=True):
    reactions_list = []

    if reactions==True :
        query = 'is_reaction'
    else:
        query = 'is_chemical'

    if query in tree.keys():
        if tree[query]:
            reactions_list.append(tree[field])
            for child in tree['children']:
                reactions_list = reactions_list + get_all_children(child, field, reactions)
    else:
        for child in tree['children']:
            reactions_list = reactions_list + get_all_children(child, field, reactions)
    return reactions_list

def get_pathway_length (result):
    return len(np.unique(get_all_children(result)))

def flatten (list_of_lists):
    return [x for ls in list_of_lists for x in ls]

def get_pathway_by_prioritizer (results, prioritizer):
    pathways = []
    for result in results['output']:
        prioritizers = []
        path_prioritizers = flatten(get_all_children(result, field = 'template_set'))
        if prioritizer in path_prioritizers:
            pathways.append(result)

    return pathways

def get_best_path (path_list, metric = 'score'):
    if len(path_list) == 0:
        return None
    if metric == 'score':
        path_scores = [np.product(get_all_children(path, 'template_score')) for path in path_list]
        best_idx = np.argmax(path_scores)
    elif metric == 'shortest':
        path_scores = [get_pathway_length(path) for path in path_list]
        best_idx = np.argmin(path_scores)
    else:
        raise ('no such metric')

    return path_list[best_idx]

def summarize_path_numbers (paths, prioritizer = None):
    num_trees = {}
    for smiles in paths:
        num_trees[smiles] = {}
        for prioritizers in paths[smiles]:
            if 'output' in paths[smiles][prioritizers].keys():
                if not prioritizer:
                    paths_ls = paths[smiles][prioritizers]['output']
                else:
                    paths_ls = get_pathway_by_prioritizer(paths[smiles][prioritizers], prioritizer)
                num_trees[smiles][prioritizers] = len(paths_ls)
            else:
                print (smiles, prioritizers, paths[smiles][prioritizers])
    return num_trees

def get_shortest_path_length (path_ls, account_for_one_pot = False, bio_templates = 'bkms'):
    best_path = get_best_path(path_ls, 'shortest')
    if best_path != None:
        shortest_path = get_all_children(best_path, 'template_set')
        if not account_for_one_pot:
            return len(shortest_path)
        else:
            path_sources = []
            for template_source in shortest_path:
                if bio_templates in template_source:
                    path_sources.append('bio')
                else:
                    path_sources.append('chem')
            #deduplicate consecutive appearances of "bio" templates
            #to simulate one-put reaction as a single step
            if len(path_sources) == 1:
                return 1

            deduplicated = [i for i, j in zip(path_sources, path_sources[1:]+[None])
                                            if i != j or i == 'chem']

            return len(deduplicated)
    else:
        print ("Could not find best path")
        return None

def summarize_shortest_paths (paths, prioritizer = None, account_for_one_pot = False):
    shortest_path_lens = {}

    for smiles in paths:
        shortest_path_lens[smiles] = {}
        for prioritizers in paths[smiles]:
            if 'output' in paths[smiles][prioritizers].keys():
                if not prioritizer:
                    paths_ls = paths[smiles][prioritizers]['output']
                else:
                    paths_ls = get_pathway_by_prioritizer(paths[smiles][prioritizers], prioritizer)
                if len(paths_ls) > 0:
                    shortest_path_lens[smiles][prioritizers] = get_shortest_path_length(paths_ls, account_for_one_pot)
            else:
                print (smiles, prioritizers, paths[smiles][prioritizers])
    return shortest_path_lens


def load_results (results_path):
    try:
        with open(results_path, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, ValueError):
        with open(results_path, 'r') as f:
            results = {}
            print(results_path)
            lines = [json.loads(entry) for entry in f.readlines()[0].replace('}{','}\t{').split('\t')]
            for line in lines:
                smiles = list(line.keys())[0]
                prioritizer = list(line[smiles].keys())[0]
                if smiles in results.keys():
                    results[smiles][prioritizer] = line[smiles][prioritizer]
                else:
                    results[smiles] = {}
                    results[smiles][prioritizer] = line[smiles][prioritizer]
            return results



def get_venn_regions(all_summary, bio_summary):
    bkms = set(bio_summary.loc['bkms',:][bio_summary.loc['bkms',:]>0].index)
    reaxys = set(all_summary.loc['reaxys',:][all_summary.loc['reaxys',:]>0].index)
    chemoenzymatic = set(bio_summary.loc['bkms,reaxys',:][bio_summary.loc['bkms,reaxys',:]>0].index)
    all_found = len(bkms.intersection(reaxys).intersection(chemoenzymatic))
    bkms_and_chemoenzymatic = len(bkms.intersection(chemoenzymatic)) - all_found
    bkms_and_reaxys = len(bkms.intersection(reaxys)) - all_found
    chemoenzymatic_and_reaxys = len(chemoenzymatic.intersection(reaxys)) - all_found
    just_bkms = len(bkms) - bkms_and_reaxys - bkms_and_chemoenzymatic - all_found
    just_reaxys = len(reaxys) - bkms_and_reaxys - chemoenzymatic_and_reaxys - all_found
    just_chemoenzymatic = len(chemoenzymatic) - bkms_and_chemoenzymatic - chemoenzymatic_and_reaxys - all_found
    return just_bkms, just_reaxys, bkms_and_reaxys, just_chemoenzymatic, bkms_and_chemoenzymatic, chemoenzymatic_and_reaxys, all_found

def get_cluster_labels(all_summary, bio_summary):
    bkms = set(bio_summary.loc['bkms',:][bio_summary.loc['bkms',:]>0].index)
    reaxys = set(all_summary.loc['reaxys',:][all_summary.loc['reaxys',:]>0].index)
    chemoenzymatic = set(bio_summary.loc['bkms,reaxys',:][bio_summary.loc['bkms,reaxys',:]>0].index)
    all_found = bkms.intersection(reaxys).intersection(chemoenzymatic)
    bkms_and_chemoenzymatic = bkms.intersection(chemoenzymatic) - all_found
    bkms_and_reaxys = bkms.intersection(reaxys) - all_found
    chemoenzymatic_and_reaxys = chemoenzymatic.intersection(reaxys) - all_found
    just_bkms = bkms - bkms_and_reaxys - bkms_and_chemoenzymatic - all_found
    just_reaxys = reaxys - bkms_and_reaxys - chemoenzymatic_and_reaxys - all_found
    just_chemoenzymatic = chemoenzymatic - bkms_and_chemoenzymatic - chemoenzymatic_and_reaxys - all_found
    none = set(all_summary.columns) - bkms - reaxys - chemoenzymatic

    label_dict = {smiles:0 for smiles in none}
    label_dict.update({smiles:1 for smiles in just_bkms})
    label_dict.update({smiles:1 for smiles in just_bkms})
    label_dict.update({smiles:2 for smiles in just_reaxys})
    label_dict.update({smiles:3 for smiles in just_chemoenzymatic})
    label_dict.update({smiles:4 for smiles in bkms_and_reaxys})
    label_dict.update({smiles:5 for smiles in bkms_and_chemoenzymatic})
    label_dict.update({smiles:6 for smiles in chemoenzymatic_and_reaxys})
    label_dict.update({smiles:7 for smiles in all_found})

    label_dict_dict = {0:'none', 1:'just bkms', 2:'just reaxys', 3: 'just chemoenzymatic',
                      4: 'bkms and reaxys', 5: 'bkms and chemoenzymatic',
                      6: 'chemoenzymatic and reaxys', 7: 'all'}

    return label_dict, label_dict_dict

def graphify_reactions_ls(reactions_ls, source_ls, numerical = {'bkms':1, 'reaxys_enzymatic':1,'reaxys':2}):

    nums = numerical #convert template names to numbers
    edges = []
    for reaction, sources in zip(reactions_ls, source_ls):
        reactants = reaction.split('>>')[0].split('.')
        product = reaction.split('>>')[1]
        numerical_source = sum([nums[source] for source in np.unique(sources)])
        edges.append((reaction, product, numerical_source))
        for reactant in reactants:
            edges.append((reactant, reaction, numerical_source))
    return edges

def process_template_source (template_source_ls, numerical = {'bkms':1, 'reaxys_enzymatic':1, 'reaxys':2}):
    nums = numerical #convert template names to numbers
    output = []
    for path in template_source_ls:
        path_output = []
        for sources in path:
            numerified = sum([nums[source] for source in np.unique(sources)])
            path_output.append(numerified)
        output.append(path_output)
    return(output)

def make_graph (target, reactions, sources, output, save_graph = True,
                numerical = {'bkms':1, 'reaxys_enzymatic':1, 'reaxys':2}, display=True, tform_dict=None):
    num_sources = process_template_source(sources, numerical = numerical) #color edges based on [numerical]
    edges = [graphify_reactions_ls(reaction_ls, sources_ls, numerical = numerical) for reaction_ls, sources_ls in zip(reactions, sources)]
    molecules = np.unique([mol for path in edges for step in path for mol in step[:2]])
    mols = [Chem.Draw.MolToImage(Chem.MolFromSmiles(smi), size=(600,600)).save('mol_images/'+smi.replace('/','(fs)').replace('\\','(bs)')+'.jpg') for smi in molecules if ('>>' not in smi) and (Chem.MolFromSmiles(smi))]
    std_target = Chem.MolToSmiles(Chem.MolFromSmiles(target)) #standardize SMILES

    edges_fl = flatten(edges)
    g = Graph()
    smiles = g.new_vp('string')
    mol_img = g.new_vp('string')
    source_e = g.new_ep('vector<float>', [''])
    vertex_size = g.new_vp('int')
    vertex_type = g.new_vp('string')
    vertex_shape = g.new_vp('string')
    vertex_fill_color = g.new_vp('string')
    vertex_frame_color = g.new_vp('string')

    if tform_dict:
        tforms_v = g.new_vp('vector<string>')

    cmap = plt.get_cmap('Set2')

    for mol in molecules:
        v = g.add_vertex()
        smiles[v] = mol

        if '>' in mol or not Chem.MolFromSmiles(mol):
            mol_img[v] = 'mol_images/black.jpg'
            vertex_size[v] = 20
            vertex_type[v] = 'reaction'
            vertex_shape[v] = 'circle'
            vertex_fill_color[v] = 'black'
            if tform_dict:
                tforms_v[v] = tform_dict[mol]
        else:
            mol_img[v] = 'mol_images/' + smiles[v].replace('/','(fs)').replace('\\','(bs)') + '.jpg'
            vertex_size[v]=90
            vertex_type[v] = 'molecule'
            vertex_shape[v] = 'square'
            vertex_fill_color[v] = 'white'
    added = []
    for edge in edges_fl:
        if (list(smiles).index(edge[0]), list(smiles).index(edge[1])) not in added:
            e = g.add_edge(list(smiles).index(edge[0]), list(smiles).index(edge[1]))
            source_e[e] = cmap(edge[2]-1)
            added.append((list(smiles).index(edge[0]), list(smiles).index(edge[1])))

    for v in g.get_vertices():
        if g.vertex(v).in_degree() == 0:
            vertex_frame_color[v] = 'green'
        elif smiles[v] == std_target or smiles[v] == target:
            vertex_frame_color[v] = 'red'
        else:
            vertex_frame_color[v] = 'black'

    g.vertex_properties['smiles'] = smiles
    g.vertex_properties['mol_img'] = mol_img
    g.vertex_properties['vertex_size'] = vertex_size
    g.vertex_properties['vertex_type'] = vertex_type
    g.vertex_properties['vertex_shape'] = vertex_shape
    g.vertex_properties['vertex_fill_color'] = vertex_fill_color
    g.vertex_properties['vertex_frame_color'] = vertex_frame_color
    g.edge_properties['source'] = source_e
    if tform_dict:
        g.vertex_properties['tform'] = tforms_v

    pos_sfdp = graph_tool.draw.sfdp_layout(g, C=1.5, K=100, gamma=0.5, p =4)
    g.vertex_properties['sfdp'] = pos_sfdp
    if display:
        graph_tool.draw.graph_draw(g, pos = pos_sfdp, vertex_surface=mol_img, output_size = (2000,2000),
                                   vertex_size=vertex_size, vertex_color=vertex_frame_color,
                                   vertex_fill_color=vertex_fill_color, vertex_shape = vertex_shape, edge_pen_width = 3,
                                   edge_color = source_e, edge_length=10, vertex_pen_width=5,
                                   nodesfirst=False, output=output)
    return g

def save_graph (outputs, smiles, prioritizers, output_prefix, buyables = 'all_buyables',
                numerical = {'bkms':1, 'reaxys_enyzmatic':1, 'reaxys':2}):
    reactions = [get_all_children(x, 'smiles') for x in outputs]
    sources = [get_all_children(x, 'template_set') for x in outputs]
    clean_smiles = smiles.replace('/','(fs)').replace('\\','(bs)')
    if len(reactions) != 0:
        g = make_graph(smiles, reactions, sources, output=output_prefix+'{}_{}_{}.png'.format(clean_smiles, buyables, prioritizers), numerical = numerical)
    else:
        return None
    return g


def show_shortest_path (result_ls):


    lengths = [get_pathway_length(result) for result in result_ls]

    min_length = np.min(lengths)
    idxs = [idx for idx, length  in enumerate(lengths) if length==min_length]

    pathways = []
    for i, idx in enumerate(idxs):
        print ('Pathway number {}'.format(i))
        reactions = [get_all_children(result) for result in result_ls][idx]
        tforms = [get_all_children(result, 'tforms') for result in result_ls][idx]
        sources = [get_all_children(result, 'template_set') for result in result_ls][idx]

        show_reaction_ls(reactions, metadata = list(zip(tforms, sources, reactions)))

        pathways.append((reactions, tforms, sources))
    return pathways


def show_paths_by_length (result_ls, query_length):


    lengths = [get_pathway_length(result) for result in result_ls]


    idxs = [idx for idx, length  in enumerate(lengths) if length==query_length]

    pathways = []
    for i, idx in enumerate(idxs):
        print ('Pathway number {}'.format(i))
        reactions = [get_all_children(result) for result in result_ls][idx]
        tforms = [get_all_children(result, 'tforms') for result in result_ls][idx]
        sources = [get_all_children(result, 'template_set') for result in result_ls][idx]

        show_reaction_ls(reactions, metadata = list(zip(tforms, sources, reactions)))

        pathways.append((reactions, tforms, sources))
    return pathways

def plot_distance_comparison_bars (shortest_paths_1, shortest_paths_2, save_path, color=True):
    improvement_dict = dict(Counter((shortest_paths_1 - shortest_paths_2).dropna()))

    thresh_number = 5

    for_pie_chart = {'≥ {} steps shorter'.format(thresh_number):0, 'Only longer':0}
    for k,v in sorted(improvement_dict.items()):
        if k <= -thresh_number:
            for_pie_chart['≥ {} steps shorter'.format(thresh_number)] += v
        elif k < -1:
            for_pie_chart['{} steps shorter'.format(str(int(-k)))] = v
        elif k == -1:
            for_pie_chart['1 step shorter'] = v
        elif k == 0:
            for_pie_chart['Equivalent length'] = v
        elif k > 0:
            for_pie_chart['Only longer'] += v

    # to get ordering right
    longer = for_pie_chart.pop('Only longer')
    for_pie_chart['Only longer'] = longer


    labels= list(for_pie_chart.keys())
    counts = list(for_pie_chart.values())

    cmap = plt.get_cmap('Greens')
    print (labels)
    print (counts)
    if color:
        colors = []
        for label in labels:
            if 'longer' in label:
                colors.append('lavenderblush')
            elif 'shorter' in label:
                steps = int(re.findall('[1-9]', label)[0])
                colors.append(cmap((steps)/(thresh_number+2)))
            else:
                colors.append('lightgray')
    else:
        colors = len(counts) * ['lightgray']
    print ("Total number of pathways compared: {}".format(np.sum(counts)))
    # plt.figure(figsize=(10,10))
    # plt.pie(np.flip(counts), labels=np.flip(labels), colors = np.flip(colors), autopct="%.1f%%", wedgeprops={'linewidth':2, 'edgecolor':'black'})
    # plt.savefig('/Users/Itai/Box Sync/Grad/Manuscripts/Chemonzymatic_planner/matplotlib_figures/boutiques_improvement_pie_chart.pdf')
    plt.figure(figsize=(2,2.5))
    plt.bar(labels, counts, color = colors, edgecolor='black')
    plt.ylabel('Number of pathways')
    plt.xticks(rotation = 90)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi =300)

    plt.show()

    
#only add reaction node to the queue if all of its reactants have been reached
#maintain a separate queue for reaction nodes that haven't yet been reached 

def metabolic_dijkstra (g, starting_nodes, account_for_one_pot=False, one_pot_template_set='bkms'):
    """
    Returns a graph where each node is labeled with the number of enzymatic 
    steps required to reach the node from a pre-defined set of starting nodes

    Implementation inspiration drawn from: 
    https://gist.github.com/aeged/db5bfda411903ecd89a3ba3cb7791a05
    
    g : (networkx.DiGraph) metabolic network
    starting_nodes : (collection) list of source nodes defined to have distance of 0
    """
    g = g.copy()
    visited = set()
    pq_1 = PriorityQueue()  
    pq_2 = PriorityQueue()  

    #initialize
    
    for node in g.nodes:
        g.nodes[node]['path_length'] = np.inf
        g.nodes[node]['visited'] = False
        g.nodes[node]['shortest_pathway'] = []
    
    for node in starting_nodes:
        g.nodes[node]['path_length'] = 0
        pq_1.put((g.nodes[node]['path_length'], node))
        visited.add(node)
        g.nodes[node]['visited'] = True
        
    while pq_1.qsize() > 0:
        curr_dist, curr_node = pq_1.get()
        # Replenish queue with reaction nodes that were not originally reached
        if pq_1.qsize()==0:
            pq_1 = pq_2

        # reaction node
        if '>>' in curr_node:
            products = list(g.successors(curr_node))
            for product in products:
                if g.nodes[product]['path_length'] > curr_dist:
                    g.nodes[product]['path_length'] = curr_dist
                    g.nodes[product]['shortest_pathway'] = g.nodes[curr_node]['shortest_pathway'] 
                if product not in visited:
                    visited.add(product)
                    pq_1.put((g.nodes[product]['path_length'], product))


        #chemical node
        else:     

            reactions = list(g.successors(curr_node))
            for reaction in reactions:
                reactants = list(g.predecessors(reaction))
                # all precursors have been visited
                if all([r in visited for r in reactants]):
#                     print (len(g.nodes[reaction]['shortest_pathway']))
                    if account_for_one_pot and len(g.nodes[curr_node]['shortest_pathway']) > 0:
                        curr_reaction_templates = g.nodes[reaction]['template_set']
                        prev_reaction = g.nodes[curr_node]['shortest_pathway'][-1]
                        prev_reaction_templates = g.nodes[prev_reaction]['template_set']
#                         print (prev_reaction_templates)
#                         print (curr_reaction_templates)
                        if one_pot_template_set in curr_reaction_templates and one_pot_template_set in prev_reaction_templates:
#                             print ('Consecutive bio steps')
                            path_length = sum([g.nodes[r]['path_length'] for r in reactants])
#                             print (path_length)
                        else:
#                             print ('Steps not consecutive')
#                             print (prev_reaction, reaction)
                            path_length = sum([g.nodes[r]['path_length'] for r in reactants]) + 1
#                             print (path_length)
                        
                    else:
#                         print ('This case')
                        path_length = sum([g.nodes[r]['path_length'] for r in reactants]) + 1 #all previous steps + 1
            
                    
                    g.nodes[reaction]['path_length'] = path_length
                    g.nodes[reaction]['shortest_pathway'] = g.nodes[curr_node]['shortest_pathway'] + [reaction]
                    pq_1.put((path_length, reaction))

                else:
                    pq_2.put((g.nodes[reaction]['path_length'], reaction))
    
    return g


####FUNCTIONS COPIED FROM ASKCOS#####
import uuid
NODE_LINK_ATTRS = {'source': 'from', 'target': 'to', 'name': 'id', 'key': 'key', 'link': 'edges'}
NIL_UUID = '00000000-0000-0000-0000-000000000000'

def clean_json(path):
    """
    Clean up json representation of a pathway. Accepts paths from either
    tree builder version.

    Note about chemical/reaction node identification:
        * For treedata format, chemical nodes have an ``is_chemical`` attribute,
          while reaction nodes have an ``is_reaction`` attribute
        * For nodelink format, all nodes have a ``type`` attribute, whose value
          is either ``chemical`` or ``reaction``
    """
    # Only fields in this dictionary are kept
    key_map = {
        'smiles': 'smiles',
        'id': 'id',
        'as_reactant': 'as_reactant',
        'as_product': 'as_product',
        'plausibility': 'plausibility',
        'purchase_price': 'ppg',
        'template_score': 'template_score',
        'terminal': 'terminal',
        'tforms': 'tforms',
        'num_examples': 'num_examples',
        'necessary_reagent': 'necessary_reagent',
        'precursor_smiles': 'precursor_smiles',
        'rms_molwt': 'rms_molwt',
        'num_rings': 'num_rings',
        'scscore': 'scscore',
        'rank': 'rank',
        'template_set': 'template_set'
    }

    if 'nodes' in path:
        # Node link format
        nodes = []
        for node in path['nodes']:
            new_node = {key_map[key]: value for key, value in node.items() if key in key_map}
            new_node['type'] = node['type']
            nodes.append(new_node)
        path['nodes'] = nodes
        output = path
    else:
        # Tree data format
        output = {}
        for key, value in path.items():
            if key in key_map:
                output[key_map[key]] = value
            elif key == 'type':
                if value == 'chemical':
                    output['is_chemical'] = True
                elif value == 'reaction':
                    output['is_reaction'] = True
            elif key == 'children':
                output['children'] = [clean_json(c) for c in value]

        if 'children' not in output:
            output['children'] = []

    return output

def nx_paths_to_json(paths, root_uuid, json_format='treedata'):
    """
    Convert list of paths from networkx graphs to json.
    """
    if json_format == 'treedata':
        # Include graph attributes at top level of resulting json
        return [{'attributes': path.graph, **clean_json(nx.tree_data(path, root_uuid))} for path in paths]
    elif json_format == 'nodelink':
        return [clean_json(nx.node_link_data(path, attrs=NODE_LINK_ATTRS)) for path in paths]
    else:
        raise ValueError('Unsupported value for json_format: {0}'.format(json_format))
        
def generate_unique_node():
    """
    Generate a unique node label using the UUID specification.

    Use UUIDv4 to generate random UUIDs instead of UUIDv1 which is used by
    ``networkx.utils.generate_unique_node``.
    """
    return str(uuid.uuid4())


def get_paths(tree, root, root_uuid=NIL_UUID, max_depth=None, max_trees=None, validate_paths=True):
    """
    Generate all paths from the root node as `nx.DiGraph` objects.

    All node attributes are copied to the output paths.

    Args:
        validate_paths (bool): require all leaves to meet terminal criteria

    Returns:
        generator of paths
    """
    import itertools

    def get_chem_paths(_node, _uuid, chem_path):
        """
        Return generator of paths with current node as the root.
        """
        if tree.out_degree(_node) == 0 or (max_depth != None and len(chem_path) >= max_depth):
            sub_path = nx.DiGraph()
            sub_path.add_node(_uuid, smiles=_node, **tree.nodes[_node])
            if (sub_path.nodes[_uuid]['terminal'] and validate_paths) or not validate_paths:
                yield sub_path
            else:
                return
        else:
            for rxn in tree.successors(_node):
                rxn_uuid = generate_unique_node()
                for sub_path in get_rxn_paths(rxn, rxn_uuid, chem_path + [_node]):
                    sub_path.add_node(_uuid, smiles=_node, **tree.nodes[_node])
                    sub_path.add_edge(_uuid, rxn_uuid)
                    yield sub_path

    def get_rxn_paths(_node, _uuid, chem_path):
        """
        Return generator of paths with current node as root.
        """
        precursors = list(tree.successors(_node))
        if len(set(precursors).intersection(set(chem_path)))>0:
            # Adding this reaction would create a cycle
            return
        c_uuid = {c: generate_unique_node() for c in precursors}
        for path_combo in itertools.product(*(get_chem_paths(c, c_uuid[c], chem_path) for c in precursors)):
            sub_path = nx.union_all(path_combo)
            sub_path.add_node(_uuid, smiles=_node, **tree.nodes[_node])
            for c in tree.successors(_node):
                sub_path.add_edge(_uuid, c_uuid[c])
            yield sub_path

    num_paths = 0
    for path in get_chem_paths(root, root_uuid, []):
        if max_trees is not None and num_paths >= max_trees:
            break
        #if validate_paths and validate_path(path) or not validate_paths:
            # Calculate depth of this path, i.e. number of reactions in longest branch
        path.graph['depth'] = [path.nodes[v]['type'] for v in nx.dag_longest_path(path)].count('reaction')
        num_paths += 1
        yield path

def prune(tree, root, max_prunes=100):
    """
    Returns a pruned networkx graph. Iteratively removes non-"terminal" leaf nodes
    and their associated parent reaction nodes.

    Args:
        tree (nx.DiGraph): full results graph from tree builder expansion
        root (str): node ID of the root node (i.e. target chemical)
        max_prunes (int): maximum number of pruning iterations. 

    Returns:
        nx.DiGraph with non-terminal leaf nodes and parent reaction nodes removed
    """
    pruned_tree = tree.copy()
    num_nodes = pruned_tree.number_of_nodes()

    for i in range(max_prunes):
        non_terminal_leaves = [
            v
            for v, d in pruned_tree.out_degree()
            if d == 0 and not pruned_tree.nodes[v]["terminal"] and v != root
        ]

        for leaf in non_terminal_leaves:
            pruned_tree.remove_nodes_from(list(pruned_tree.predecessors(leaf)))
            pruned_tree.remove_node(leaf)

        if pruned_tree.number_of_nodes() == num_nodes:
            break
        else:
            num_nodes = pruned_tree.number_of_nodes()

    # If pruning resulted in a disconnected graph, remove nodes not connected to the root subgraph
    if not nx.is_weakly_connected(pruned_tree):
        for c in nx.weakly_connected_components(pruned_tree):
            if root in c:
                pruned_tree.remove_nodes_from([n for n in pruned_tree if n not in c])
                break

    return pruned_tree

