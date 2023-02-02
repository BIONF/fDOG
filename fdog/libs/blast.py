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

import os
import sys
from Bio.Blast.Applications import NcbiblastpCommandline
import xml.etree.ElementTree as ET
import subprocess


##### FUNCTIONS RELATED TO BLAST #####

def do_blastsearch(
        query, blast_db, evalBlast = 0.00001, lowComplexityFilter = False):
    """ Perform blastp search for a query fasta file
    Return an XML string contains blast result
    """
    filter = 'no'
    if lowComplexityFilter == True:
        filter = 'yes'
    try:
        blastp_cline = NcbiblastpCommandline(
            query = query, db = blast_db, evalue = evalBlast, seg = filter,
            max_target_seqs = 10, outfmt = 5)
        stdout, stderr = blastp_cline()
        return(stdout)
    except:
        sys.exit(
            'ERROR: Error running blastp search for %s against %s'
            % (query, blast_db))


def parse_blast_xml(blast_output):
    """ Parse Blast XML output from a string variable
    Return a dictionary containing query ID, query length, together with a list
    of hits and their bit score, evalue, align len
    """
    blast_dict = {}
    root = ET.fromstring(blast_output)
    blast_dict['query'] = root[8][0][2].text
    blast_dict['query_len'] = root[8][0][3].text
    blast_dict['hits'] = {}
    for type_tag in root.findall(
            'BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
        value = type_tag.findall('*')
        hit_id = 'NA'
        for i in type_tag.findall('*'):
            if i.tag == 'Hit_def':
                if not i.text in blast_dict['hits']:
                    hit_id = i.text
                    blast_dict['hits'][hit_id] = {}
            if i.tag == 'Hit_hsps':
                if hit_id in blast_dict['hits']:
                    blast_dict['hits'][hit_id]['bit_score'] = i[0][1].text
                    blast_dict['hits'][hit_id]['evalue'] = i[0][3].text
                    blast_dict['hits'][hit_id]['align_len'] = i[0][13].text
    return(blast_dict)


def make_blastdb(args):
    """ Make blastDB in coreTaxa_dir
    for fdog.addTaxon, fdog.addTaxa and fdog.checkData
    """
    (specName, specFile, outPath, coreTaxa_dir, searchTaxa_dir, silent) = args
    if not coreTaxa_dir:
        coreTaxa_dir = '%s/coreTaxa_dir' % outPath
    if not searchTaxa_dir:
        searchTaxa_dir = '%s/searchTaxa_dir' % outPath
    blastCmd = 'makeblastdb -dbtype prot -in %s -out %s/%s/%s' % (specFile, coreTaxa_dir, specName, specName)
    if silent == True:
        blastCmd = blastCmd + '> /dev/null 2>&1'
    try:
        subprocess.call([blastCmd], shell = True)
    except:
        sys.exit('Problem with running %s' % blastCmd)
    fileInGenome = "%s/%s/%s.fa" % (searchTaxa_dir, specName, specName)
    fileInBlast = "%s/%s/%s.fa" % (coreTaxa_dir, specName, specName)
    if not os.path.exists(fileInBlast):
        os.symlink(fileInGenome, fileInBlast)
