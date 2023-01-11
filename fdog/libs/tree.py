# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This file is part of fDOG tool https://github.com/BIONF/fDOG
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: tran@bio.uni-frankfurt.de
#
#######################################################################


from ete3 import NCBITaxa

import fdog.libs.zzz as general_fn


##### FUNCTIONS RELATED TO TAXONOMY TREE #####

def get_rank_index(lineage, rank_name, ncbi):
    """ Get ID and index in the species lineage for a given rank
    Return {rank_id:rank_index}
    """
    ranks = ncbi.get_rank(lineage)
    rank_id = list(general_fn.matching_elements(ranks, rank_name).keys())[0]
    rank_index = len(ranks) - lineage.index(rank_id) - 1
    return({rank_id:rank_index})


def get_rank_range(lineage, minDist, maxDist, ncbi):
    """ Get rank ID and its index in a given species lineage
    for a pair of min and max rank. See get_rank_index()
    Return a list of 2 dictionary for min and max rank as
    [{min_rank_id:min_rank_index}, {max_rank_id:max_rank_index}]
    """
    return(
        get_rank_index(lineage, minDist, ncbi),
        get_rank_index(lineage, maxDist, ncbi))


def check_taxon_group(group_id, tax_id, ncbi):
    """ Check if a taxon (tax_id) belongs to a taxonomy group (group_id)"""
    lineage = ncbi.get_lineage(tax_id)
    if group_id in lineage:
        return(True)
    return(False)


def get_ancestor(id1, id2, ncbi):
    """ Get common ancestor ID and rank for 2 taxon IDs
    Return dictionary {ancestor_id: ancestor_rank}
    """
    tree = ncbi.get_topology([id1, id2], intermediate_nodes = False)
    ancestor = tree.get_common_ancestor(id1, id2).name
    return(ncbi.get_rank([ancestor]))


def check_common_ancestor(ref_id, ancestor, minDist, maxDist, ncbi):
    """ Check if ancestor ID lies within the range between min and max rank
    of reference species
    Return 1 if true
    """
    ref_lineage = ncbi.get_lineage(ref_id)
    (min_ref, max_ref) = get_rank_range(ref_lineage, minDist, maxDist, ncbi)
    ancestor_index = len(ref_lineage) - ref_lineage.index(ancestor) - 1
    if list(min_ref.values())[0] <= ancestor_index <= list(max_ref.values())[0]:
        return(1)
    return(0)


def remove_clade(tree, node_id):
    """ Remove a clade from a tree """
    removed_clade = tree.search_nodes(name = str(node_id))[0]
    removed_node = removed_clade.detach()
    return(tree)


def get_leaves_dict(spec_lineage, tree, min_index, max_index):
    """ Given a tree and a lineage string of a species
    Return a dictionary where keys are the internal nodes defined by the
    ranks between min rank (e.g. genus, specified by min_index in the species
    lineage) and max rank (e.g. phylum). Values are all leaves in the tree
    that belong to the corresponding internal node (rank)
    """
    node_dict = {}
    already_added = []
    spec_lineage.reverse()
    for i in range(len(spec_lineage)):
        if i >= min_index and i <= max_index:
            curr_node = spec_lineage[i]
            node = tree.search_nodes(name = str(curr_node))
            if len(node) > 0:
                for leaf in node:
                    node_dict[spec_lineage[i]] = []
                    for t in leaf.traverse():
                        if t.is_leaf():
                            if not t.name in already_added:
                                already_added.append(t.name)
                                node_dict[spec_lineage[i]].append(t.name)
    return(general_fn.remove_dup_in_dict(node_dict))


def get_tax_name(taxId):
    """ Get taxonomy name for a given taxon ID """
    ncbi = NCBITaxa()
    try:
        ncbiName = ncbi.get_taxid_translator([taxId])[int(taxId)]
        ncbiName = re.sub('[^a-zA-Z1-9\s]+', '', ncbiName)
        taxName = ncbiName.split()
        name = taxName[0][:3].upper()+taxName[1][:2].upper()
    except:
        name = "UNK" + taxId
    return(name)
