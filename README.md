# Tree_topology_compare  

## Purpose

Compare tree topologies and score shared branching. When comparing tree topology or branching most available tools score "concensus" between trees of identical leaves or taxons. Therefore, in addition to consensus, it can also score branching between trees of partially shared taxons, e.g., Jack-knife Monopyly Indexing (JMI), to determine the impact of taxon sampling in tree topology.

## Process

1. Collect trees (tree file format supported by Dendropy, e.g., Newick or Nexus) from given file paths as an argument.  
2. Reroot trees.  
3. Find and compare topology of shared taxons between two trees.  
4. Determine a shared branching between a reference tree and other trees, using two scoring types by default.

   * H1: shared branching between trees of identical taxons, e.g., concensus.  
   * H2: shared branching between trees of partially shared taxons, e.g., Jack-Knife Monophyly Indexing (JMI).  

## Requirements

* Dendropy (4.5.2)  
* Python3  
Lower version of libraries or packages may also support this script but not tested.  

## Usage

python3 [this_scipt_path] [arguments] [tree_file_path(s)] > save_file_path  

### Input

Path to tree files, either in Newick or Nexus format, and also accept a file contain multiple trees.  

### Output

Standard output (print to screen) branch scored tree in Nexus format. Unfortunately, Newick format cannot fully or properly show metadata. In linux or terminal system, use '>' to store output. 

### Arguments

* -h, show help()
* -r [path], path to the reference tree to be scored and printed to screen
* -n, normalize scores (count / a number of other trees)
* -R [str], reroot using one or more taxons, which is necessary, unless all input trees are rerooted to an idential taxon.
  Use ',' (comma) as a delimiter to input more than one taxons ("taxon1,taxon2" -> taxon1, taxon2).
  Given more than one taxons will be rerooted by their most recent common ancestor node.
* -m [str], dicard 'H1' or 'H2' score names and manually use one given score name
* -f [str], tree (file) format or schema. Dendropy support Newick or Nexus tree schema
