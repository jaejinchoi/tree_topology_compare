### Tree_topology_compare  

### Purpose
Compare tree topologies and score shared branching. When comparing tree topology or branching most available tools score "concensus" between trees of identical operational taxon units (OTUs). In addition to determine consensus, this script can also score branching of partially shared OTUs, e.g., Jack-knife Monopyly Indexing (JMI), to determine the impact of OTU sampling in tree topology.

### Process
1. Collect trees (tree file format supported by Dendropy, e.g., Newick or Nexus) from given file paths and as an argument.
2. Reroot trees
3. Keep and compare shared operation taxon units (OTUs) between two trees for comparison.
4. Score branchings by comparing topology between trees, a reference tree vs. other trees, using two scoring types by default.
  * H1: shared branching between trees of identical OTUs, e.g., concensus.
  * H2: shared branching between trees of partially shared OTUs, e.g., Jack-Knife Monophyly Indexing (JMI).

  
### Requirements
* Dendropy (4.5.2)
* Python3
Lower version of libraries or packages may also support this script but not tested. 


### Input
Path to tree files, either in Newick or Nexus format, and also accept a file contain multiple trees.

### Output
Standard output (print to screen) branch scored tree in Nexus format. Unfortunately, Newick format cannot fully or properly show metadata. In linux or terminal system, use '>' to store output. 

### Usage
python3 [this_scipt_path] [arguments] [tree_file_path(s)] > save_file_path
  
## Arguments
* -h, show help()
* -r [path], path to the reference tree to be scored and printed to screen
* -n, normalize scores (count / a number of other trees)
* -R [str], reroot using one or more OTUs, which is necessary, unless all input trees are preroot to an idential OTU.
  Use ',' (comma) as a delimiter to input more than one OTUs ("OTU1,OTU2" -> OTU1, OTU2).
  Given more than one OTUs will be rerooted by their most common recent ancestor node.
* -m [str], dicard 'H1' or 'H2' score names and manually use one given score name
* -f [str], tree (file) format or schema. Dendropy support Newick or Nexus tree schema

