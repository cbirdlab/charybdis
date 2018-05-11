#!/usr/bin/python
import os
import math
import re
import numpy as nm
import pandas as pd
import xml.etree.ElementTree
import pickle
import SAP
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import urllib

def similarityScore(seq1, seq2):
# Source: SAP codebase, 'UtilityFunctions.py'
    assert (len(seq1) == len(seq2))

    seq1LeftBoundary = len(re.search("^(-*)", seq1).groups()[0])
    seq1RightBoundary = len(seq1) - len(re.search("(-*)$", seq1).groups()[0])
    seq2LeftBoundary = len(re.search("^(-*)", seq2).groups()[0])
    seq2RightBoundary = len(seq2) - len(re.search("(-*)$", seq2).groups()[0])
    leftBoundary = max(seq1LeftBoundary, seq2LeftBoundary)
    rightboundary = min(seq1RightBoundary, seq2RightBoundary)

    seq1Trunc = seq1[leftBoundary:rightboundary]
    seq2Trunc = seq2[leftBoundary:rightboundary]

    length = len(seq1Trunc)

    extraGaps = 0
    ident = 0
    for i in range(0, length):
        if seq1Trunc[i] == "-" and seq2Trunc[i] == "-":
            extraGaps += 1
            continue
        if seq1Trunc[i] == seq2Trunc[i]:
            ident += 1
    score = ident/float(length - extraGaps)
    return score


def getPathFromSEQID(SEQID, target):
    # target: The type of file requested,
    # such as a treestatscache file
    return 0

def getBranchSummaryRecurse(node):
    children = node.keys()
    numChildren = len(children)
    if numChildren == 0:
        return {'count':0, 'depth':1, 'leaves':[node], 'branchLeafCount':0} 
    summary = {}
    summary['count'] = numChildren - 1
    summary['branchLeafCount'] = 0
    summary['depth'] = 1
    summary['leaves'] = []
    depths = []
    for child in children:
        ret = getBranchSummaryRecurse(node[child])
        summary['leaves'] = summary['leaves'] + ret['leaves']
        summary['count'] = summary['count'] + ret['count']
        depths.append(ret['depth'])
        summary['branchLeafCount'] = ret['branchLeafCount'] + summary['branchLeafCount']

    summary['depth'] = summary['depth'] + max(depths)
    if (all(x == 1 for x in depths) == True):
        summary['branchLeafCount'] = numChildren - 1
    return summary
def getBranchSummary(tree):
    # dTree is a tree stored as a dictionary
    root = tree[tree.keys()[0]]
    summary = {}
    summary['count'] = 1
    summary['branchLeafCount'] = 0
    summary['depth'] = -1
    summary['leaves'] = []
    ret = getBranchSummaryRecurse(root) 
    summary['count'] = summary['count'] + ret['count']
    summary['depth'] = summary['depth'] + ret['depth'] 
    summary['leaves'] = ret['leaves']
    summary['branchLeafCount'] = ret['branchLeafCount']
    return summary

def getAssignmentsFile(SAPdir):
    assignmentsPath = SAPdir + "assignments.csv"
    print (assignmentsPath)
    #assignmentsPath = SAPdir + "assignments-2.csv"
    if (os.path.exists(assignmentsPath) == False):
        print ("File not found: ", assignmentsPath)
        return None
    assignments = pd.read_csv(assignmentsPath, header=None)
    assignments.columns = ["File", "Cutoff", "Detail", "ID", 
            "Phylum", "Class", "Order", "Family", "Genus", "Species",
            "NumHomologues", "MinFreqHomologue",  "MinTaxonProb"]
    return assignments

def getTaxonProbabilitiesFile(SAPdir):
    taxonProbsPath = SAPdir + "taxon_probabilities.csv"
    if (os.path.exists(taxonProbsPath) == False):
        print ("File not found: ", taxonProbsPath)
        return None
    taxonProbs = pd.read_csv(taxonProbsPath, header=None)
    taxonProbs.columns = ["File", "ID", "BestRank", "BestTaxon", "PosteriorProb"]
    return taxonProbs

def getTreeFile(SAPdir, queryFile, SEQID):
    treePath = SAPdir + "treestatscache/" + queryFile + "_" + SEQID.replace(":", "_").replace("-", "_").replace(" ", "_").replace("|", "_") + ".pickle"
    if (os.path.exists(treePath) == False):
        print ("File not found: ", treePath)
        return None
    treeFile = open(treePath, 'rb')
    tree = pickle.load(treeFile)
    treeFile.close()
    tree._removePlaceHolders()
    return tree

def getBlastCacheFile(SAPdir, queryFile, SEQID):
    # MUST be dynamic/paramterized
    blastCachePath = SAPdir + "blastcache/" + queryFile + "_" + SEQID.replace(":", "_").replace("-", "_").replace(" ", "_").replace("|", "_") + ".200_10.0" + ".xml"
    if (os.path.exists(blastCachePath) == False):
        print ("File not found: ", blastCachePath)
        return None
    blastCache = xml.etree.ElementTree.parse(blastCachePath).getroot()
    return blastCache

def getHomologCacheFile(SAPdir, queryFile, SEQID):
    homologCachePath = SAPdir + "homologcache/" + queryFile + "_" + SEQID.replace(":", "_").replace("-", "_").replace(" ", "_").replace("|", "_") + ".pickle"
    if (os.path.exists(homologCachePath) == False):
        print ("File not found: ", homologCachePath)
        return None
    homologCacheFile = open(homologCachePath, "rb")
    homologCache = pickle.load(homologCacheFile)
    homologCacheFile.close()
    return homologCache

def getHomologFastaFile(SAPdir, queryFile, SEQID):
    homologFastaPath = SAPdir + "homologcache/" + queryFile + "_" + SEQID.replace(":", "_").replace("-", "_").replace(" ", "_").replace("|", "_") + ".fasta"
    if (os.path.exists(homologFastaPath) == False):
        print ("File not found: ", homologFastaPath)
        return None
    homologFasta = SeqIO.to_dict(SeqIO.parse(homologFastaPath, "fasta"))
    return homologFasta

def getClonesFile(SAPdir, queryFile, SEQID):
    clonesPath = SAPdir + "html/clones/" + queryFile + "_"  + SEQID.replace(":", "_").replace("-", "_").replace(" ", "_").replace("|", "_") + ".html"
    if (os.path.exists(clonesPath) == False):
        print ("File not found: ", clonesPath)
        return None
    clones = urllib.urlopen(clonesPath).read()
    return clones

def ensure_dir(file_path):
    #Source: https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def makeHomologueFastaPath(SAPdir, queryFile, SEQID):
    ensure_dir(SAPdir + "evaluation/")
    homologueFastaPath = SAPdir + "evaluation/" + "homologuesToAnalyze_" + queryFile + "_" + SEQID.replace(":", "_").replace("-", "_").replace(" ", "_").replace("|", "_") + ".fasta"
    return homologueFastaPath

def makeBlastResultsPath(SAPdir, queryFile, SEQID):
    ensure_dir(SAPdir + "evaluation/")
    blastResultsPath = SAPdir + "evaluation/" + "blast_" + queryFile + "_" + SEQID.replace(":", "_").replace("-", "_").replace(" ", "_").replace("|", "_") + ".xml"
    return blastResultsPath

def makeOutputPath(SAPdir, outfile):
    outputPath = SAPdir + outfile + ".csv"
    return outputPath

def most_common(lst):
    # Source: https://stackoverflow.com/a/1518632
    return max(set(lst), key=lst.count)


def main():
    # Parameters
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-d", "--directory", dest="DIR",
            help = "SAP results directory (root)")
    parser.add_option("-c", "--cutoff", dest="CUT", 
            help = "evaluate results at this cutoff")
    parser.add_option("-o", "--outfile", dest="OUT",
            help = "Output file to hold attributes CSV")
    (options, args) = parser.parse_args ()
    
    # SHOULD BE PARAMETERS
    #SAPdir = "/media/Wapuilani/evan/repo/SAP_formatter_clean/Simons8/"
    SAPdir = options.DIR
    cutoff = float (options.CUT)
    outfile = options.OUT


    if (os.path.isdir(SAPdir) == False):
        print ("Directory", SAPdir, "does not exist")
        exit(1)
    if (SAPdir[len(SAPdir) - 1] != '/'):
        SAPdir = SAPdir + "/"
    
    # Get assignments file
    assignments = getAssignmentsFile(SAPdir)

    assignments = assignments.loc[assignments['Cutoff'] == cutoff]

    # Get taxon probabilities file
    taxonProbs = getTaxonProbabilitiesFile(SAPdir)
    taxonProbsBest = taxonProbs.drop_duplicates(subset="ID")

    # Initialize table to hold scenario classification for each query sequence
    classified = taxonProbsBest['ID'].copy().to_frame()
    classified['SpeciesLevel'] = False
    classified['SpeciesCount'] = 0
    classified['BranchCount'] = 0
    classified['BranchLeafCount'] = 0
    classified['BranchInteriorCount'] = 0
    classified['TreeDepth'] = 0
    classified['SpeciesProbMax'] = 0
    classified['SpeciesProbMin'] = 0
    classified['HomologuesToAnalyzeCount'] = 0
    classified['HomologuesSuspiciousRecordsCount'] = 0
    classified['HasHomologuesSuspiciousRecords'] = False
    classified['HomologueCount'] = 0
    #classified['MinFreqHomologueProb'] = 0
    #classified['MinFreqHomologueProbWarning'] = False
    #classified['MinTaxonProb'] = 0
    #classified['MinTaxonProbWarning'] = False
    classified['DatabaseNotExhausted'] = False


    # Investigate and classify each sequence assignment
    # Note that the number of sequences may by less than the number of query sequences,
    # since some queries cannot be assigned given the parameters and database
    for index, row in taxonProbsBest.iterrows():

        row['ID'] = re.sub (r" $", "", row['ID'])
        print row['ID']
        #break

        print (assignments.loc[assignments['ID'] == row['ID'], 'Species'])
        # Is there a species-level identification (posterior prob >= cutoff)
        try:
            math.isnan(assignments.loc[assignments['ID'] == row['ID'], 'Species'])
            classified.loc[classified['ID'] == row['ID'], 'SpeciesLevel'] = False
        except:
            classified.loc[classified['ID'] == row['ID'], 'SpeciesLevel'] = True


        # Get number of branches and tree depth
        tree = getTreeFile (SAPdir, row['File'], row['ID'])
        branchSummary = getBranchSummary(tree)
        classified.loc[classified['ID'] == row['ID'], 'BranchCount'] = branchSummary['count']
        classified.loc[classified['ID'] == row['ID'], 'TreeDepth'] = branchSummary['depth']

        classified.loc[classified['ID'] == row['ID'], 'BranchLeafCount'] = branchSummary['branchLeafCount']
        classified.loc[classified['ID'] == row['ID'], 'BranchInteriorCount'] = branchSummary['count'] - branchSummary['branchLeafCount']
        # Get number of species in tree
        speciesCount = len(branchSummary['leaves'])
        classified.loc[classified['ID'] == row['ID'], 'SpeciesCount'] = speciesCount

        # Get max and min species posterior prob
        speciesProbs = tree.getLevelProbs(level="species")[-speciesCount:]
        print (speciesProbs)
        if (len (speciesProbs) > 0):
            species, probs = zip(*speciesProbs)
            speciesProbsMax = max(probs)
            classified.loc[classified['ID'] == row['ID'], 'SpeciesProbMax'] = speciesProbsMax
            speciesProbsMin = min(probs)
            classified.loc[classified['ID'] == row['ID'], 'SpeciesProbMin'] = speciesProbsMin
        else:
            classified.loc[classified['ID'] == row['ID'], 'SpeciesProbMax'] = float('NaN')
            classified.loc[classified['ID'] == row['ID'], 'SpeciesProbMin'] = float('NaN')


        # Get HomologueCount, MinFreqHomologueProb, MinTaxonProb
        
        classified.loc[classified['ID'] == row['ID'], 'HomologueCount'] = len(assignments.loc[assignments['ID'] == row['ID'], 'NumHomologues'])
        #classified.loc[classified['ID'] == row['ID'], 'MinFreqHomologueProb'] = float(assignments.loc[assignments['ID'] == row['ID'], 'MinFreqHomologue'])
        #if classified.loc[classified['ID'] == row['ID'], 'MinFreqHomologueProb'].any() > 0.0001:
        #    classified.loc[classified['ID'] == row['ID'], 'MinFreqHomologueProbWarning'] = True
        #classified.loc[classified['ID'] == row['ID'], 'MinTaxonProb'] = float (assignments.loc[assignments['ID'] == row['ID'], 'MinTaxonProb'])
        #if classified.loc[classified['ID'] == row['ID'], 'MinTaxonProb'].any()  > 0.0001:
        #    classified.loc[classified['ID'] == row['ID'], 'MinTaxonProbWarning'] = True

        # BLAST each species in tree

        MislabelledHomologueCount = 0

        blastCache = getBlastCacheFile(SAPdir, row['File'], row['ID'])
        homologCache = getHomologCacheFile(SAPdir, row['File'], row['ID'])
        homologFasta = getHomologFastaFile(SAPdir, row['File'], row['ID'])
        homologues = homologCache.homologues
        clones = getClonesFile(SAPdir, row['File'], row['ID'])

        # Check for warning about database not being exhausted
        warning = re.search("database was not exhausted", clones)
        if warning is not None:
            classified.loc[classified['ID'] == row['ID'], 'DatabaseNotExhausted'] = True


        querySequence = homologCache.queryFasta.sequence
        
        homologuesToAnalyze = []
        homologuesToAnalyzeFasta = []
        homologuesToAnalyzeFastaFile = makeHomologueFastaPath(SAPdir, row['File'], row['ID'])
        homologuesToAnalyzeBlastResultsFile = makeBlastResultsPath(SAPdir, row['File'], row['ID'])

        for h in homologues:
            hFastaHeader = h + "_" + h
            hFastaEntry = homologFasta[hFastaHeader]
            cloneLine = re.search(hFastaHeader + r"<\/a>:<\/td><td>[.0-9]*<\/td>", clones)
            hIdentityScore = re.search(r"[0-9]*\.[0-9]*", cloneLine.group(0)).group(0)
            if (float(hIdentityScore) > 0.9):
                homologuesToAnalyze.append(dict({'ID':h, 'FastaEntry':hFastaEntry, 'IdenityScore':hIdentityScore}))
                homologuesToAnalyzeFasta.append(hFastaEntry)

        classified.loc[classified['ID'] == row['ID'], 'HomologuesToAnalyzeCount'] = len(homologuesToAnalyze)


        # If the number of homologues to analyze is greater than 1, Blast to find mislabeled database entries
        if (len(homologuesToAnalyze) > 1):

            # Write selected homologues to Fasta in order to Blast them
            SeqIO.write(homologuesToAnalyzeFasta, homologuesToAnalyzeFastaFile, "fasta")

            # Blast them
            blastx_cline = NcbiblastnCommandline(query = homologuesToAnalyzeFastaFile, 
                    db="/media/Wapuilani/Databases/FILTER_BLASTDB/nr_mito.SAP.fix.fasta", outfmt=5, 
                    out=homologuesToAnalyzeBlastResultsFile, perc_identity=90)
            blastx_cline()

            result = open(homologuesToAnalyzeBlastResultsFile, "r")
            records = NCBIXML.parse(result)
            hasRecords = True

            # Investigate Blast results for each selected homologue
            item = next(records)
            suspiciousHomologueDatabaseRecords = 0
            while (hasRecords == True):
                genusList = []
                queryGenus = ""
                count = 0
                for alignment in item.alignments:
                    try:
                        genus = re.search(r"genus[^,]*,", alignment.title).group(0).replace('genus: ', '').replace(',', '')
                    except:
                        genus = ""
                    genusList.append(genus)
                    if count == 0:
                        queryGenus = genus
                    count = count + 1

                consensusGenus = most_common(genusList)
                if (consensusGenus != queryGenus):
                    suspiciousHomologueDatabaseRecords = suspiciousHomologueDatabaseRecords + 1

                item = next(records, False)
                if (item == False):
                    hasRecords = False

            classified.loc[classified['ID'] == row['ID'], 'HomologuesSuspiciousRecordsCount'] = suspiciousHomologueDatabaseRecords
            if suspiciousHomologueDatabaseRecords > 0:
                classified.loc[classified['ID'] == row['ID'], 'HasHomologuesSuspiciousRecords'] = True

        #print classified.loc[classified['ID'] == row['ID']]

        #if (len(homologuesToAnalyze) > 1):
        #   exit(0)

    outputPath = makeOutputPath(SAPdir, outfile)
    classified.to_csv(outputPath, index = False)

    return 0

if __name__ == "__main__":
    main()
