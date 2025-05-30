# Functions for programs call of the basics MSA / phylogenetic tools

import os
import sys
import shutil
import subprocess
import codecs
import configparser
from Bio import SeqIO
from ete3 import Tree
from pathlib import Path

import dot_fastaBlocs
import agraph_fastaBlocs

# Get the way to call the programm from the config file

config_file = f"{os.path.abspath(os.path.dirname(__file__))}/config.txt"
config = configparser.ConfigParser()
config.readfp(codecs.open(config_file, "r", "utf-8-sig"))



def msa(fasta_file) -> tuple:
    """
    Launch a multi sequences alignement
    """
    out_fn = Path(f"{fasta_file.parents[0]}/muscle_{fasta_file.stem}.fasta").resolve()
    if config['ENV']['MUSCLE']:
        env = config['ENV']['MUSCLE'].split(' ') 
    else: 
        env = []
    cmd = [
            config['PROGRAMS']['MUSCLE'],
            "-in", f"{fasta_file}",
            "-out", f"{out_fn}",
            "-quiet"
            ]
    cmd = ' '.join(env + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def all_msa(directory) -> tuple:
    """
    Alignement for all the fasta_file of a given directory
    """
    process_list = []
    for filename in Path(directory).iterdir():
        filename.resolve()
        if filename.suffix == ".fasta":
            process, out_fn = msa(filename)
            process_list.append(process)
    return process_list

def trimal(msa_fasta, gt="0.9", cons="05") -> tuple:
    """
    Colonns selection of a MSA using trimal
    """
    out_fn = Path(f"{msa_fasta.parents[0]}/trimal_{msa_fasta.stem}.fasta").resolve()
    if config['ENV']['TRIMAL']:
        env = config['ENV']['TRIMAL'].split(' ') 
    else: 
        env = []
    cmd = [
            config['PROGRAMS']['TRIMAL'],
            "-in", f"{msa_fasta}",
            "-out", f"{out_fn}",
            "-gt", gt,
            "-cons", cons
            ]
    cmd = ' '.join(env + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def phylo(msa_fasta, d="aa", extra_args_phyml=None) -> tuple:
    """
    Phylogenetic inference of a given msa_file
    """
    node = []
    with open(msa_fasta, "r") as f_file:
        for line in f_file:
            if line.startswith(">"):
                node.append(line.replace(">","").replace("\n",""))
    phylip_file = Path(f"{msa_fasta.parents[0]}/{msa_fasta.stem}.phylip").resolve()
    records = SeqIO.parse(msa_fasta, "fasta")
    SeqIO.write(records, phylip_file, "phylip-relaxed")
    out_fn = Path(f"{phylip_file}_phyml_tree.txt").resolve()
    if config['ENV']['PHYML']:
        env = config['ENV']['PHYML'].split(' ') 
    else: 
        env = []
    cmd = [
            config['PROGRAMS']['PHYML'],
            "-i", f"{phylip_file}",
            "-d", d,
            "--no_memory_check"
            ]
    # Add extra arguments if provided
    if extra_args_phyml:
        cmd += extra_args_phyml
    cmd = ' '.join(env + cmd)
    if len(node) > 2:
        process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    elif len(node) == 2:
        with open(out_fn, "w+") as t_file:
            t_file.write(f"({node[0]},{node[1]});")
        process = subprocess.Popen("echo -n", shell=True, stdout=open(os.devnull, 'wb'))
    elif len(node) == 1:
        process = subprocess.Popen("echo -n", shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def all_phylo(directory, extra_args_phyml=None) -> tuple:
    """
    Phylogenetic inference for all the msa_file of a given directory
    """
    process_list = []
    for filename in Path(directory).iterdir():
        filename.resolve()
        if filename.suffix == ".fasta":
            process, out_fn = phylo(filename, extra_args_phyml=extra_args_phyml)
            process_list.append(process)
    return process_list

def all_msa_phylo(directory) -> tuple:
    """
    Alignement for all the fasta_file of a given directory
    """
    process_list = []
    for filename in Path(directory).iterdir():
        filename.resolve()
        if filename.suffix == ".fasta":
            process, out_fn = msa(filename)
            process_list.append(process)
    [p.wait() for p in process_list]
    process_list = []
    for filename in Path(directory).iterdir():
        filename.resolve()
        if filename.suffix == ".fasta" and filename.stem.startswith("muscle_"):
            process, out_fn = phylo(filename)
            process_list.append(process)
    return process_list

def compute_branch_length(msa_fasta, treefix_tree, d="aa", o="lr") -> tuple:
    """
    Use maximum likelihood compute branch length of a fixed topology tree
    (based on the initial MSA)
    """
    node = []
    with open(msa_fasta, "r") as f_file:
        for line in f_file:
            if line.startswith(">"):
                node.append(line.replace(">","").replace("\n",""))
    phylip_file = Path(f"{msa_fasta.parents[0]}/{msa_fasta.stem}.phylip").resolve()
    records = SeqIO.parse(msa_fasta, "fasta")
    SeqIO.write(records, phylip_file, "phylip-relaxed")
    out_fn = Path(f"{phylip_file}_phyml_tree.txt").resolve()
    if config['ENV']['PHYML']:
        env = config['ENV']['PHYML'].split(' ')
    else:
        env = []
    cmd = [
            config['PROGRAMS']['PHYML'],
            "-i", f"{phylip_file}",
            "-d", d,
            "-u", f"{treefix_tree}",
            "-o", o,
            "--no_memory_check"
            ]
    cmd = ' '.join(env + cmd)
    if len(node) > 2:
        process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    elif len(node) == 2:
        with open(out_fn, "w+") as t_file:
            t_file.write(f"({node[0]},{node[1]});")
        process = subprocess.Popen("echo -n", shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def segmentation(fasta_file, q="2", m="1", M="20", t="10", m_iter=0, extra_args_paloma=None) -> tuple:
    """
    Module segmentation using paloma-2 (partial local multiple alignment)
    Output file being the .agraph file (yaml format)
    """
    if config['ENV']['PALOMA-2']:
        env = config['ENV']['PALOMA'].split(' ') 
    else: 
        env = []

    # Init list for not already used paloma arg
    e_extra_args_paloma = []
    # Parse extra args
    if extra_args_paloma:
        extra_split = extra_args_paloma
        i = 0
        while i < len(extra_split):
            arg = extra_split[i]
            val = extra_split[i + 1] if i + 1 < len(extra_split) else ""
            # Override if already defined
            if '-t' or '--thres' in arg:
                t = int(val)
            elif '-q' or '--quorum' in arg:
                q = int(val)
            elif '-m' or '--min-size' in arg:
                m = int(val)
            elif '-M' or '--max-size' in arg:
                M = int(val)
            else:
                e_extra_args_paloma.append(f"{arg} {val}")
            i += 2  # move to next pair

    out_fn = Path(f"{fasta_file.stem}_t{t}m{m}M{M}.oplma").resolve()
    # If not first iteration, just return the name of the oplma that is already computed
    if int(m_iter) > 0:
        process = subprocess.Popen([sys.executable, "-c", "exit(0)"], stdout=open(os.devnull, 'wb'))
        return process, out_fn

    cmd = [
            config['PROGRAMS']['PALOMA-2'],
            "-i", f"{fasta_file}",
            "-q", f'{q}',
            "-m", f'{m}',
            "-M", f'{M}',
            "-t", f'{t}',
            "--oplma",
            ]
    # Add extra arguments if provided
    if len(e_extra_args_paloma) > 0:
        cmd += e_extra_args_paloma
    cmd = ' '.join(env + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def modules_fasta(ms_output) -> str:
    """
    Create a directory containing modules fasta file, for a given paloma output file
    """
    module_directory = Path(f"{ms_output.parents[0]}/modules_{ms_output.stem}").resolve()
    if not os.path.exists(module_directory):
            os.makedirs(module_directory)
    if ms_output.suffix == ".agraph":
        agraph_fastaBlocs.make_module_directory(ms_output, module_directory)
    elif ms_output.suffix == ".oplma":
        # Build dot file from old plma
        output_dot = Path(f"{ms_output.parents[0]}/{ms_output.stem}.dot").resolve()
        if config['ENV']['PLMA2DOT']:
            env = config['ENV']['PLMA2DOT'].split(' ') 
        else: 
            env = []
        cmd = [
                config['PROGRAMS']['PLMA2DOT'],
                "-i", f"{ms_output}",
                "-o", f"{output_dot}",
                ]
        cmd = ' '.join(env + cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
        process.wait()
        thres = 5
        dot_fastaBlocs.make_module_directory(output_dot, thres, module_directory)
    elif ms_output.suffix == ".dot":
        # Directly from a dot
        thres = 5
        dot_fastaBlocs.make_module_directory(ms_output, thres, module_directory)
    return module_directory

def treefix(msa_fasta, phylo_tree, species_tree, A=".fasta", model="PROTGAMMAJTT", V="0", niter="100", extra_args_treefix=None, extra_args_raxml=None) -> tuple:
    """
    Correction of a given tree, using a guide tree
    """
    treefix_dir = Path(f"{phylo_tree.parents[0]}/{phylo_tree.stem}_treefix_dir").resolve()
    if not os.path.exists(treefix_dir):
        os.makedirs(treefix_dir)
    shutil.copy(msa_fasta, Path(f"{treefix_dir}/{phylo_tree.stem}.fasta").resolve())
    shutil.copy(phylo_tree, Path(f"{treefix_dir}/{phylo_tree.stem}.tree").resolve())
    smap_file = make_smap(treefix_dir, species_tree)
    tree_to_fix_path = f"{treefix_dir.parents[0]}/{treefix_dir.name}_treeToFixPath.txt"
    os.system(f"ls -d -1 {treefix_dir}/*tree > {tree_to_fix_path}")
    out_fn = Path(f"{treefix_dir}/{phylo_tree.stem}.treefix.tree").resolve()
    if config['ENV']['TREEFIX']:
        env = config['ENV']['TREEFIX'].split(' ') 
    else: 
        env = []
    
    # Set default RAxML args and allow override of niter
    raxml_args = []
    dtl_args = []
    # Parse extra args
    if extra_args_raxml:
        extra_split = extra_args_raxml
        i = 0
        while i < len(extra_split):
            arg = extra_split[i]
            val = extra_split[i + 1] if i + 1 < len(extra_split) else ""
            if "-m" in arg:
                model = val  # override model
            else:
                raxml_args.append(f"{arg} {val}")
            i += 2  # move to next pair
    raxml_args.append(f"-m {model}")
    if extra_args_treefix:
        extra_split = extra_args_treefix
        i = 0
        while i < len(extra_split):
            arg = extra_split[i]
            val = extra_split[i + 1] if i + 1 < len(extra_split) else ""
            if "--niter" in arg:
                niter = val  # override niter
            else:
                dtl_args.append(f"{arg} {val}")
            i += 2  # move to next pair

    cmd = [
            config['PROGRAMS']['TREEFIX'],
            "-i", f"{tree_to_fix_path}",
            "-s", f"{species_tree}",
            "-S", f"{smap_file}",
            "-A", str(A),
            "-V", str(V),
            "--niter", f"{niter}",
            "-e", f"'{' '.join(raxml_args)}'",
            "-E", f"'{' '.join(dtl_args)}'"
            ]
    cmd = ' '.join(env + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def make_smap(directory, specie_tree) -> str:
    """
    Read all the trees in the directory, make them binaries, then write all their mappings
    in the smap file
    """
    outName = Path(f"{directory.parents[0]}/{directory.stem}.smap").resolve()
    smapFile = open(outName, "w+")
    sp_tree = Tree(f"{specie_tree}")
    for file in Path(directory).iterdir():
        if file.suffix == ".tree":
            treeFile = open(file, "r")
            for line in treeFile:
                try:
                    tree =  Tree(line, format = 0)
                except:
                    print(f"a file is going to be deleted, check it ! {file}")
                    os.remove(file)
                    break
                tree.resolve_polytomy(recursive=True)
                tree.write(format = 9, outfile = file)
                for node in tree.traverse("postorder"):
                    node.name = node.name.replace("  ", "")
                    # If there is a | in the name, it's a module tree
                    if "|" in node.name:
                        # So in this case, the specie is the gene
                        specie = "_".join(node.name.split("_")[1:-1])
                    # If not, it's a gene tree
                    else:
                        specie = node.name.split("_")[-1]
                    # The specie / gene must corrspond to a node in the specie / gene tree    
                    for sp_node in sp_tree.traverse("postorder"):
                        if specie in sp_node.name and len(specie) > 2:
                            smapFile.write(f"{node.name}\t{sp_node.name}\n")
            treeFile.close()    
    smapFile.close()
    return outName

def pastml(gene_tree, pastml_csv, sep=',', extra_args=None) -> tuple:
    """
    Run pastML
    """
    out_fn = Path(f"{pastml_csv.parents[0]}/{pastml_csv.stem}_combined_ancestral_states.tab").resolve()
    if config['ENV']['PASTML']:
        env = config['ENV']['PASTML'].split(' ')
    else:
        env = []
    cmd = [
            config['PROGRAMS']['PASTML'],
            "-t", f"{gene_tree}",
            "-d", f"{pastml_csv}",
            "-s", sep,
            "-o", f"{out_fn}",
            #"--forced_joint"
            ]
    # Add extra arguments if provided
    if extra_args:
        cmd += extra_args
    cmd = ' '.join(env + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def seadog_md(species_tree, gene_tree, path_modules_tree, m_iter=0, extra_args=None) -> tuple:
    """
    Make a DGS phylogenetic reconciliation, using SEADOG-MD
    """
    out_fn = Path(f"{m_iter+1}_seadogMD_{gene_tree.stem.split('_')[-1].split('.')[0]}.output").resolve()
    gene_tree_directory = Path(f"{species_tree.parents[0]}/gene_tree_{gene_tree.stem.split('.')[0]}").resolve()
    if not os.path.exists(gene_tree_directory):
        os.makedirs(gene_tree_directory)
    shutil.copy(gene_tree, Path(f"{gene_tree_directory}/{gene_tree.stem.split('_')[-1].split('.')[0]}.tree").resolve())
    if config['ENV']['SEADOG-MD']:
        env = config['ENV']['SEADOG-MD'].split(' ') 
    else: 
        env = []
    cmd = [
            config['PROGRAMS']['SEADOG-MD'],
            "-s", f"{species_tree}",
            "-g", f"{gene_tree_directory}/",
            "-d", f"{path_modules_tree}",
            "-o", f"{out_fn}",
            ]
    # Add extra arguments if provided
    if extra_args:
        cmd += extra_args
    cmd = ' '.join(env + cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn

def known_domains(multi_fasta) -> tuple:
    """
    Search for known domains and motif on our proteins sequences of interests
    """
    get_domains_program = f"{os.path.abspath(os.path.dirname(__file__))}/get_domains.py" 
    out_fn = Path(f"{multi_fasta.parents[0]}/domains_{multi_fasta.stem}.csv").resolve()
    cmd = [ 
            "python3",
            get_domains_program,
            "--fasta", f"{multi_fasta}",
            "-o", f"{out_fn}"
        ]
    cmd = ' '.join(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=open(os.devnull, 'wb'))
    return process, out_fn


