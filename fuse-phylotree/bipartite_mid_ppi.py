# Build itol file for a expression file

import argparse
from pathlib import Path
import random
import networkx as nx
from networkx.algorithms import bipartite
import plotly.graph_objects as go
import seaborn as sns
from ete3 import Tree

#==============================================================================
# Read files
#==============================================================================

def read_csv(csv_fn) -> dict:
    """ 
    Read csv file 
    format : gene,ppi_gain,ppi_lost,mod_gain,mod_lost
    return dict 
    {gene : ppi_gain, ppi_lost, mod_gain, mod_lost}
    """
    dict_gene_association = {}
    with open(csv_fn) as csv_file:
        for line in csv_file:
            if "gene" in line:
                continue
            else:
                splited_line = line.replace("\n","").split(",")
                gene = splited_line[0]
                ppi_gain_list = [i.split("_")[1] for i in splited_line[1].split("|") if len(i) > 1]
                ppi_lost_list = [i.split("_")[1] for i in splited_line[2].split("|") if len(i) > 1]
                mod_gain_list = [i for i in splited_line[3].split("|") if len(i) > 1]
                mod_lost_list = [i for i in splited_line[4].split("|") if len(i) > 1]
                dict_gene_association[gene] = [ppi_gain_list, ppi_lost_list, mod_gain_list, mod_lost_list]
    return dict_gene_association

#==============================================================================
# Filters and operation
#==============================================================================

def co_appeared(dict_gene_association) -> dict:
    """
    Filter the gene association dict
    Keep only gene/association where at least 1 module and 1 ppi are gained
    Return subdict with same format
    """
    return {
        gene : association 
        for gene, association in dict_gene_association.items()
        if len(association[0]) > 0 and len(association[2]) > 0
        }

#==============================================================================
# Graph operation
#==============================================================================

def bipartite_graph(dict_gene_association) -> object:
    """
    Build bipartite association graph
    representing at on side : modules, at the other side ppi
    an edge represents a gene node where module and ppico-appeared (association)
    """
    G = nx.Graph()
    for gene, association in dict_gene_association.items():
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        G.add_nodes_from(mod_gain_list, bipartite=0)
        G.add_nodes_from([gene], bipartite=1)
        G.add_nodes_from(ppi_gain_list, bipartite=2)
        for ppi in ppi_gain_list:
            G.add_edges_from([(gene, ppi)], gene=gene)   
        for mod in mod_gain_list:
            G.add_edges_from([(gene, mod)], gene=gene)    
    assert bipartite.is_bipartite(G)
    return G

def full_graph(dict_gene_association, gene_tree) -> object:
    """
    Full association graph
    nodes : gene, modules or PPI
    edge : gained at this gene OR a gene is a descendend of another
    """
    G = nx.Graph()
    for gene, association in dict_gene_association.items():
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        G.add_nodes_from(mod_gain_list, bipartite="module")
        G.add_nodes_from([gene], bipartite="gene")
        G.add_nodes_from(ppi_gain_list, bipartite="PPI")
        for ppi in ppi_gain_list:
            G.add_edges_from([(gene, ppi)], gene=gene)   
        for mod in mod_gain_list:
            G.add_edges_from([(gene, mod)], gene=gene)    
        #gene_node = gene_tree&gene
    # If gene is a direct descendant (children) of another, add a edge
    for gene_node in gene_tree.traverse():
        if gene_node.name not in dict_gene_association:
            G.add_nodes_from([gene_node.name], bipartite="gene")
        for child_node in gene_node.children:
            G.add_edges_from([(gene_node.name, child_node.name)], phylo_rel="Children")
    return G

def random_hexaColor() -> str:
    """
    Generate an string with a random hexadecimal color code
    """
    col = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    return col

def sankey_diag(dict_gene_association) -> object:
    """
    Build Sankey diagram using plotly package
    """
    k = 0
    dict_elem_nb = {}
    label, color = [], []
    palette = sns.color_palette("deep", len(dict_gene_association))
    palette = [f"rgb({c[0]},{c[1]},{c[2]})" for c in palette]
    i = 0
    for gene, association in dict_gene_association.items():
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        if gene not in dict_elem_nb:
            dict_elem_nb[gene] = k
            label.append(gene)
            color.append(palette[i])
            i += 1
            k += 1
        for ppi in ppi_gain_list:
            if ppi not in dict_elem_nb:
                dict_elem_nb[ppi] = k
                label.append(ppi)
                color.append("rgb(100,149,237)")
                k += 1
        for mod in mod_gain_list:
            if mod not in dict_elem_nb:
                dict_elem_nb[mod] = k
                label.append(mod)
                color.append("rgb(50,205,50)")
                k += 1
    node = dict(
        pad = 1,
        thickness = 100,
        label = label,
        color = color
        )
    source, target, value, color = [], [], [], []
    i = 0
    for gene, association in dict_gene_association.items():
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        col = palette[i]
        i += 1
        for mod in mod_gain_list:
            source.append(dict_elem_nb[mod])
            label.append(mod)
            target.append(dict_elem_nb[gene])
            value.append(1)
            color.append(col)
        for ppi in ppi_gain_list:
            source.append(dict_elem_nb[gene])
            target.append(dict_elem_nb[ppi])
            value.append(1)
            color.append(col)
    link = dict(
        source = source, 
        target = target, 
        value = value, 
        color = color
        )
    data=go.Sankey(
        node=node, 
        link=link
        )
    fig = go.Figure(data)
    fig.add_annotation(dict(font=dict(color="rgb(50,205,50)",size=12), x=0, y=1.06, showarrow=False, text='<b>Modules</b>'))
    fig.add_annotation(dict(font=dict(color="black",size=12), x=0.5, y=1.06, showarrow=False, text='<b>Gene where module(s) and PPI(s) co-appeared</b>'))
    fig.add_annotation(dict(font=dict(color="rgb(100,149,237)",size=12), x=1, y=1.06, showarrow=False, text='<b>PPI</b>'))
    fig.update_layout(font=dict(size=5, color='black'))
    fig.show()
    fig.write_html("sankey.html")
    
def simple_sankey_diag(dict_gene_association) -> object:
    """
    Build Sankey diagram using plotly package
    """
    k = 0
    dict_elem_nb = {}
    label, color = [], []
    palette = sns.color_palette("deep", len(dict_gene_association))
    palette = [f"rgb({c[0]},{c[1]},{c[2]})" for c in palette]
    i = 0
    for gene, association in dict_gene_association.items():
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        if gene not in dict_elem_nb:
            dict_elem_nb[gene] = k
            label.append(gene)
            color.append(palette[i])
            i += 1
            k += 1
        for ppi in ppi_gain_list:
            if ppi not in dict_elem_nb:
                dict_elem_nb[ppi] = k
                label.append(ppi)
                color.append("rgb(100,149,237)")
                k += 1
        for mod in mod_gain_list:
            if mod not in dict_elem_nb:
                dict_elem_nb[mod] = k
                label.append(mod)
                color.append("rgb(50,205,50)")
                k += 1
    node = dict(
        pad = 1,
        thickness = 100,
        label = label,
        color = color
        )
    source, target, value, color = [], [], [], []
    i = 0
    for gene, association in dict_gene_association.items():
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        col = palette[i]
        i += 1
        for mod in mod_gain_list:
            source.append(dict_elem_nb[mod])
            label.append(mod)
            for ppi in ppi_gain_list:
                target.append(dict_elem_nb[ppi])
                value.append(1)
                color.append(col)
    link = dict(
        source = source, 
        target = target, 
        value = value, 
        color = color
        )
    data=go.Sankey(
        node=node, 
        link=link
        )
    fig = go.Figure(data)
    fig.add_annotation(dict(font=dict(color="rgb(50,205,50)",size=12), x=0, y=1.06, showarrow=False, text='<b>Modules</b>'))
    fig.add_annotation(dict(font=dict(color="black",size=12), x=0.5, y=1.06, showarrow=False, text='<b>Gene where module(s) and PPI(s) co-appeared</b>'))
    fig.add_annotation(dict(font=dict(color="rgb(100,149,237)",size=12), x=1, y=1.06, showarrow=False, text='<b>PPI</b>'))
    fig.update_layout(font=dict(size=5, color='black'))
    fig.show()
    fig.write_html("sankey.html")
    
def parallel_plot(dict_gene_association):
    """
    Build parallel association plot
    """
    dict_module_nb = {}
    dict_gene_nb = {}
    dict_ppi_nb = {}
    dict_gene_color = {}
    col_scale = []
    m, g, p = 0, 0, 0
    palette = sns.color_palette("deep", len(dict_gene_association))
    palette = [f"rgb({c[0]},{c[1]},{c[2]})" for c in palette]
    for gene, association in dict_gene_association.items():
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        if gene not in dict_gene_nb:
            dict_gene_nb[gene] = g
            dict_gene_color[gene] = palette[g]
            g += 1
        for ppi in ppi_gain_list:
            if ppi not in dict_ppi_nb:
                dict_ppi_nb[ppi] = p
                p += 1
        for mod in mod_gain_list:
            if mod not in dict_module_nb:
                dict_module_nb[mod] = m
                m += 1   
    df = {"module" : [], "gene" : [], "ppi" : []}
    g = 0
    for gene, association in dict_gene_association.items():
        col_scale.append([(1/len(dict_gene_nb))*g, palette[g]])
        g += 1
        ppi_gain_list = association[0]
        mod_gain_list = association[2]
        for ppi in ppi_gain_list:
            for mod in mod_gain_list:
                df["module"].append(dict_module_nb[mod])
                df["gene"].append(dict_gene_nb[gene])
                df["ppi"].append(dict_ppi_nb[ppi])
    fig = go.Figure(data=
        go.Parcoords(
            line = dict(color = df['gene'],
                        colorscale = col_scale
                        ),
            dimensions = list([
                dict(range = [0,len(dict_module_nb)],
                     ticktext = [i for i in dict_module_nb.keys()],
                     tickvals = [i for i in range(0, len(dict_module_nb))],
                     label = "Modules", values = df['module']),
                dict(range = [0,len(dict_gene_nb)],
                     tickvals = [i for i in range(0, len(dict_gene_nb))],
                     ticktext = [i for i in dict_gene_nb.keys()],
                     label = 'Genes', values = df['gene']),
                dict(range = [0,len(dict_ppi_nb)],
                     tickvals = [i for i in range(0, len(dict_ppi_nb))],
                     ticktext = [i for i in dict_ppi_nb.keys()],
                     label = 'PPIs', values = df['ppi'])
                        ])
                    )
                )
    fig.update_layout(font=dict(size=10, color='black'))
    fig.write_html("parallel.html")
    fig.show()


#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file",
                        help = "csv file containing the gene where modules and PPI co-appeared (format : gene,ppi_gain,ppi_lost,mod_gain,mod_lost, eg; G71_69_70,Q5XPI4_RNF123,,,B335|B1315)",
                        type=str)
    parser.add_argument("--gene_tree",
                        help = "newick tree file containing the evolutionary relation shared by the gene in input csv_file",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    # Get arguments
    args = parser()
    csv_fn = Path(args.csv_file).resolve()
    dict_gene_association = read_csv(csv_fn)
    print(f"{len(dict_gene_association)} module(s)-PPI(s) associations")
    if args.gene_tree:
        gene_tree = Tree(str(args.gene_tree), format=1)
        G_full = full_graph(dict_gene_association, gene_tree)
        nx.write_gexf(G_full, f"assoc_completeGraph.gexf")
        nx.write_graphml(G_full, f"assoc_completeGraph.graphml")
    dict_gene_association = co_appeared(dict_gene_association)
    print(f"{len(dict_gene_association)} co-apparitions associations")
    print(f"{len([gene for gene in dict_gene_association if gene.startswith('G')])} co-apparitions associations at ancestral gene node")
    

    for gene, assoc in dict_gene_association.items():
        if gene.startswith("G"):
            print(f"{gene.split('_')[0]} & DESC & \\tiny{{{' '.join(assoc[0])}}} & \\tiny{{{' '.join([m.replace('B', 'M') for m in assoc[2]])}}} \\\\")
            print("\midrule")
    
    print(caca)
    
    ppi_list = []
    for gene, coapp in dict_gene_association.items():
        if gene.startswith("G"):
            for ppi in coapp[0]:
                if ppi not in ppi_list:
                    ppi_list.append(ppi)
    print(f"\t{len(ppi_list)} PPI associated with module signatures")
    mod_list = []
    for gene, coapp in dict_gene_association.items():
        if gene.startswith("G"):
            for mod in coapp[2]:
                if mod not in mod_list:
                    mod_list.append(mod)
    print(f"\t{len(mod_list)} modules present in these module signatures")
    dict_gene_association = {gene : assoc for gene, assoc in dict_gene_association.items() if gene.startswith("G")}
    G = bipartite_graph(dict_gene_association)
    nx.write_gexf(G, f"assoc_bipartite.gexf")
    nx.write_graphml(G, f"assoc_mod_gene_ppi.graphml")
    sankey_diag(dict_gene_association)
    simple_sankey_diag(dict_gene_association)
    parallel_plot(dict_gene_association)
    
if __name__ == "__main__":
    main()