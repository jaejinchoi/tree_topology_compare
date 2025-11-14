#python3
import os, sys, getopt
import dendropy
#import pandas as pd

sys.setrecursionlimit(550000) #for dendropy rerooting process (recursion search limit)

class dendro_comp:
        
    def collect_taxon_list(self, tree_ob=dendropy.Tree(), exclude_taxon_list=[]): #return a list of taxons
        taxon_list = [n_node.taxon.label for n_node in tree_ob.leaf_node_iter()] #obtain subtree taxons (labels); tree.leaf_iter() deprecated
        taxon_list = list(set(taxon_list) - set(exclude_taxon_list)) #exclude unshared taxons
        taxon_list.sort()

        return taxon_list


    def subtree_component_tile_dict(self, tree_ob=dendropy.Tree(), exclude_taxon_list=[]): #return a list of subtree strings
        subtree_tile_dict={} #{subtree_tile:subtree.as_string()}

        for inner_node in tree_ob.internal_nodes():
            sub_tree = tree_ob.extract_tree(
                node_filter_fn = None #clone entire structure
                #, suppress_unifurcation=False
            ) #shallow copy (whole tree)?

            inner_node.parent_node = None #set the parent_node to None
            sub_tree.seed_node = inner_node #sub_tree.seed_node is a reference to the original tree
            #sub_tree.seed_node.parent_node = None #new seed_node with the parent_node is None
            sub_tree.update_bipartitions() #necessary per function process

            subtree_taxon_list = dendro_comp().collect_taxon_list(tree_ob=sub_tree, exclude_taxon_list=exclude_taxon_list)
            subtree_taxon_tile = ":".join(subtree_taxon_list).strip()

            #avoid end-resulting duplicated tiles and empty tiles
            if (subtree_taxon_tile!="" and subtree_taxon_tile not in subtree_tile_dict): 
                #subtree_tile_dict[subtree_taxon_tile]=sub_tree.as_newick_string()
                subtree_tile_dict[subtree_taxon_tile]=sub_tree.as_string("newick")

        return subtree_tile_dict


    def reroot_tree(self, tree=dendropy.Tree(), reroot_clade_list=[], rooting_method=""):

        if (reroot_clade_list!=[]):

            #tree = dendropy.Tree.get_from_string(newick_str, "newick")
            tree.encode_bipartitions()
            #tree.is_rooted = False #tree is unrooted state, root state place "[&R]_" at the front
            tree.is_rooted = True

            #https://dendropy.readthedocs.io/en/v3.12.1/tutorial/treemanips.html
            #reroot options: edge, node, and outgroup (single)
            if (len(reroot_clade_list)>1):

                mrca = tree.mrca(taxon_labels=reroot_clade_list) #mrca = most recent common ancestor of multiple outgroups
                
                if (rooting_method=="node"):
                    tree.reroot_at_node(mrca, update_bipartitions=True)
                
                elif (rooting_method=="edge"):
                    half_edge = mrca.edge_length / 2.0
                    tree.reroot_at_edge(mrca.edge, length1=half_edge, length2=half_edge, update_bipartitions=True, suppress_unifurcations=False)
                    #tree.reroot_at_edge(mrca.edge, update_bipartitions=True)

                else:
                    print(("unsupported rooting method using outgroups of ", reroot_clade_list))

            else: #pull one outgroup, then halved its edge
                #print(tree.taxon_namespace)
                outgroup_node = tree.find_node_with_taxon_label(reroot_clade_list[0])

                if (rooting_method=="node"):
                    tree.to_outgroup_position(outgroup_node, update_bipartitions=True)
                
                elif (rooting_method=="edge"):
                    half_edge = outgroup_node.edge_length / 2.0
                    tree.reroot_at_edge(outgroup_node.edge, length1=half_edge, length2=half_edge, update_bipartitions=True, suppress_unifurcations=False)
                
                elif (rooting_method=="midpoint"):
                    tree.reroot_at_midpoint(update_bipartitions=True) #find the longest edge and then split (not halved)

                else:
                    print(("unsupported rooting method using an outgroup of ", reroot_clade_list))

            tree.is_rooted = True

 
    def annotate_node(
        self
        , ref_tree_ob=dendropy.Tree()
        , comp_tile_dict={}
        , exclude_taxon_list=[]
        , annotation_label=""
        , annotation_score_base=0
        ):

        ## https://dendropy.org/primer/working_with_metadata_annotations.html

        for inner_node in ref_tree_ob.internal_nodes():
            #child_taxon_list = [n_leaf.taxon.label for n_leaf in inner_node.leaf_node_iter()] #raise error when inner_node==None
            child_taxon_list = [n_leaf.taxon.label for n_leaf in inner_node.leaf_iter()]

            child_taxon_list = list(set(child_taxon_list) - set(exclude_taxon_list))
            child_taxon_list.sort()

            child_taxon_tile = ":".join(child_taxon_list)

            if (inner_node.annotations.find(name=annotation_label)==None):
                inner_node.annotations[annotation_label]._value=annotation_score_base #begining from 0 or 1(self)

            if (child_taxon_tile in comp_tile_dict): #for any subtrees (inner nodes)
                #print(child_taxon_tile)
                inner_node.annotations[annotation_label]._value+=1 #additive
                #print(inner_node.annotations[annotation_label]._value)
                

def show_help():
    print('[options][load path]')
    print('-h, show_help')
    print('-v, show_description')
    print('-----input a reference tree to be scored---------')
    print('-r [path], path of reference tree to be scored')
    #print('-t, print trend of confidence score, and skip tree print')
    print('-n, normalize (count / a number of subtrees) annotation score, e.g., consensus_cnt, jmi_cnt')
    
    print("# configurations")
    print('to use without -r, add reference tree at the end of the trees file(input)')
    print('-R [str], reroot using one or more items; string delimited using a comma (e.g., item1,item2 = [item1, item2])')
    print('\tUse "," (comma) as a delimiter to input more than one taxon ("Taxon1,Taxon2" -> [Taxon1, Taxon2]')
    print('\tGiven more than one taxon will be rerooted by their most common recent ancestor node')          
    print('-m [str], a name for annotation label in an output tree')
    print('-f [str], input tree format that are supported by Dendropy (default: newick)')
    print("\tNewick")
    print("\tNexus")

    print("")
    print('# Default branch scoring types')
    print('\tH1: shared branching between trees of identical taxons (e.g., consensus)')
    print('\tH2: shared branching between trees of partially shared taxons (e.g., Jack-Knife Monophyly Indexing; JMI)')

    print("")
    print("-o [path], save a copy of (Nexus) reference tree with branch scores as a file, or standard output as default")
    print("\tDefault: %s" % (os.path.abspath(os.getcwd())+"/compared_tree.result"))
    
    sys.exit()
    

def show_version():
    print("Compare tree branching between trees of")
    print("\t1. sharing exactly same taxons or leaves (e.g., consensus)")
    print("\t2. partially sharing same taxons or leaves (e.g., Jack-knife Monophyly Indexing or JMI)")
    print("")
    print("Code by JaeJin Choi. March 6, 2023-")
    sys.exit()


if __name__=="__main__":

    #process
    '''
    1. collect trees from argument or file path
    2. keep taxons that are shared between two trees for comparison 
    3. score branching to a given reference tree [-r]
    '''
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hvr:naR:L:f:o:')

    except:
        show_help()

    ref_tree_path=''
    #save_folder_path=''
    normalize_score_flag=False

    reroot_clade_list=[] #can be one or more than two items to set as an outgroup (reroot)
    tree_format="newick"
    pre_annotation_label=""
    save_path = f"{os.getcwd()}/compared_tree.result"

    for opt, arg in opts:

        if (opt=='-h'):
            show_help()

        elif (opt=='-v'):
            show_version()

        elif (opt=='-r'):
            ref_tree_path = os.path.abspath(arg)

        elif (opt=='-n'):
            normalize_score_flag=True

        elif (opt=='-f'):
            tree_format=str(arg)

        elif (opt=='-L'): #
            pre_annotation_label=str(arg)

        elif (opt=='-R'): #rerooting
            reroot_clade_list = str(arg).split(',')

        elif (opt=='-o'): #save file path
            save_path = os.path.abspath(arg)

        else:
            print("undefined argument or flag: %s" % (opt))
            sys.exit()
            #show_help()

    exclude_taxon_list=[]
    tree_list=[]
    group_tile_dict={}

    ref_tree_ob = dendropy.Tree()
    ref_taxon_list = []

    comp_tree_ob = dendropy.Tree()
    comp_taxon_list = []

    rooting_method="node"
    annotation_label="" #consensus_cnt or jmi_cnt
    annotation_score_base=0       

    if (len(args)>0):

        if (ref_tree_path==''):
            ref_tree_path = os.path.abspath(args.pop())
            print("Reference unspecified; thus, using the last tree path in args: %s" % (ref_tree_path))

        try:
            ref_tree_ob = dendropy.Tree.get_from_path(ref_tree_path, tree_format)

        except:
            print("Unable to read a reference tree: %s" % (ref_tree_path))
            sys.exit(0)

        ## remove any duplicates of ref_tree_path in comp_tree_paths
        load_path_list = [os.path.abspath(n_path) for n_path in args]
        
        if (ref_tree_path in load_path_list):
            load_path_list.remove(ref_tree_path)

        dendro_comp().reroot_tree(
            tree=ref_tree_ob
            , reroot_clade_list=reroot_clade_list
            , rooting_method=rooting_method
            )

        ref_taxon_list = dendro_comp().collect_taxon_list(
            tree_ob=ref_tree_ob
            , exclude_taxon_list=[]
            )

        total_comp_tree_cnt=0

        for comp_tree_path in load_path_list:
            ## https://dendropy.org/primer/treecollections.html
            comp_treelist_ob = dendropy.TreeList.get(path=comp_tree_path, schema=tree_format) #read multiple trees
            #print(len(comp_treelist_ob))

            total_comp_tree_cnt+=len(comp_treelist_ob)

            for comp_tree_ob in comp_treelist_ob:
            
                dendro_comp().reroot_tree(
                    tree=comp_tree_ob
                    , reroot_clade_list=reroot_clade_list
                    , rooting_method=rooting_method
                    )

                comp_taxon_list = dendro_comp().collect_taxon_list(
                    tree_ob=comp_tree_ob
                    , exclude_taxon_list=[]
                    )

                ## exclude_taxon_list; not shared by reference and compared trees; union(A,B) - intersection(A,B)
                exclude_taxon_list = list(set(ref_taxon_list).symmetric_difference(set(comp_taxon_list)))

                print(f"# Shared taxon between REF and COMP trees | intersection : union = {len(set(ref_taxon_list) & set(comp_taxon_list))} : {len(set(ref_taxon_list) | set(comp_taxon_list))}")

                ## a name of annotation/metadata
                if (pre_annotation_label==""): #use predefined annotation label if not given

                    if (exclude_taxon_list==[]):
                        annotation_label="H1" #compare topology of trees of all shared taxons (e.g., consensus)

                    else:
                        annotation_label="H2" #compare topology between trees of partially shared taxons (e.g., Jack-knife monophyly indexing)

                else:
                    annotation_label = pre_annotation_label

                comp_tile_dict = dendro_comp().subtree_component_tile_dict(
                    tree_ob=comp_tree_ob
                    , exclude_taxon_list=exclude_taxon_list
                    )

                dendro_comp().annotate_node(
                    ref_tree_ob=ref_tree_ob
                    , comp_tile_dict=comp_tile_dict
                    , exclude_taxon_list=exclude_taxon_list
                    , annotation_label=annotation_label #eval(), locals()[], difficult to apply
                    , annotation_score_base=annotation_score_base
                    )

        print(f"\nTotal # of subtrees compared N = {total_comp_tree_cnt}, Normalize: {normalize_score_flag}")

        if (normalize_score_flag==True):
            
            for inner_node in ref_tree_ob.internal_nodes():
                if (inner_node.annotations.find(name=annotation_label)!=None):
                    #inner_node.annotations[annotation_label]._value=float(inner_node.annotations[annotation_label]._value) / float(len(comp_treelist_ob)) #begining from 0 or 1(self)
                    inner_node.annotations[annotation_label]._value=float(inner_node.annotations[annotation_label]._value) / total_comp_tree_cnt #begining from 0 or 1(self)
                    #print(inner_node.annotations[annotation_label]._value)
                        
                        
        if (save_path!=""):
            with open(save_path, 'w') as write_f:
                write_f.write(ref_tree_ob.as_string(schema="nexus"))
        
        else:
            #print((ref_tree_ob.as_string(schema=tree_format))) #cannot fully/property show (multiple) metadata using newick
            print((ref_tree_ob.as_string(schema="nexus")))
