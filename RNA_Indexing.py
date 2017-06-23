import os, re, sys, StringIO, time

#### Calling pipeline modules
import myvep
import mybasic, mysetting
sys.path.append('/data1/home/jsgene/JK1/NGS/mutation/')
import vep_batch
from GS_overload import *
from GS_download import GS_Downloads
from GS_mapping import GS_Mapping
from GS_variants import GS_calling_variants
import mymysql as msql

#### for process control
import commands as CMDs
import subprocess as sp
import Queue as qu
from multiprocessing import Process, Queue, JoinableQueue, Pool, TimeoutError
import threading as thread
import stopwatch as stw

#### for Data-mining
import numpy as np
import pandas as pd
import scipy

try:
    import GEOparse as GEO
except RuntimeError, rer:
    import GEOparse as GEO

try:
    import matplotlib as mpl
except RuntimeError, ere:
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    # plt.switch_backend('Agg')
    import matplotlib.image as mplimage
    import matplotlib.transforms as mtrans
    from matplotlib import transforms
    from matplotlib.colors import ListedColormap
from sklearn import datasets
from sklearn.neighbors import NearestCentroid
from sklearn.neighbors import KNeighborsClassifier

# Biopython
import pysam, vcf
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from reportlab.lib import colors
from reportlab.lib.units import cm, inch
from reportlab.lib.colors import magenta, red, blue, green, black, cyan, pink
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from Bio.Graphics import GenomeDiagram, BasicChromosome
from Bio import SeqIO, motifs
import weblogolib as wLogo
import corebio.seq as SEQUENCE

# output Files
import pickle, json, logging, gzip
from PyPDF2 import PdfFileMerger, PdfFileReader, PdfFileWriter
#import pylab as pl

from GS_Indexing import *
import pdfkit, wkhtmltopdf

from glob import glob
base = GS_Logging()
WORKDESK = base.Desktop
SAVE_TO = base.Save_file_path
SEQ_DIR = '/data1/Sequence/'
GENBANK = SEQ_DIR + 'Genbank/'
GENE_DICT = GENBANK + "gene_dictionary.pkl"

class Setting_for_datamining(object):
    def __init__(self):
        self.eql8 = '/EQL8/pipeline'

    def Initialise_Dataset_by_Group(self, projectN):
        groupName = self.eql8 + '/' + projectN + '_rsq2expr'
        dataFormat = groupName + "/*/*.rpkm"
        return glob(dataFormat)

    def Labeling_for_Supervised_Learning(self, trainingSet, Label):
        # setting labels
        LabelList = []
        for i in range(len(trainingSet)):
            LabelList.append(Label)
        return LabelList

    def Labeling_for_Subtype(self, subtype_Label):
        LabelList = []
        for i in range(len(subtype_Label)):
            if subtype_Label[i] == 'P':
                LabelList.append(0)
            elif subtype_Label[i] == 'N':
                LabelList.append(1)
            elif subtype_Label[i] == 'M':
                LabelList.append(2)
            elif subtype_Label[i] == 'C':
                LabelList.append(3)
            elif subtype_Label[i] == 'U':
                LabelList.append(4)
            else:
                LabelList.append(4)
        return LabelList



    def RPKM_data(self, paths):
        DF = pd.DataFrame.from_csv(paths, sep='\t', header=None)
        DF = DF.reset_index()
        DF = pd.DataFrame(DF)
        SR = pd.Series(list(DF[2]), index=list(DF[0]))
        return SR

    def Merge_DataFrame(self, fileList, dbname):

        DF = []
        ID_List = []
        for afile in fileList:

            ## From Files
            if '.rpkm' in afile:
                se = self.RPKM_data(afile)
                ID = afile.split('/')[-1]   ## e.g : '/EQL8/pipeline/SGI20160831_rsq2expr/IRCR_BT16_1031_T01_RSq/IRCR_BT16_1031_T01_RSq.rpkm'
                ID = ID.split('.rpkm')[0]   ## e.g : 'IRCR_BT16_1031_T01_RSq.rpkm'
            ## From mySQL
            else:
                dataset_from_DB = cohort_Dataset()
                #se = dataset_from_DB.get_data_from_SQL('rpkm_gene_expr', afile, dbname)
                se, ID = dataset_from_DB.get_data_from_SQL('rpkm_gene_expr', afile, dbname)
            DF.append(se)
            ID_List.append(ID)
            newDF = pd.concat(DF, axis=1)
            newDF.columns = newDF.iloc[0]

            try:
                newDF = newDF.drop(["geneName"])
            except ValueError, verr:
                continue
        print "Creating DataFrame with dimensions ", newDF.shape

        # DF_16 = newDF.iloc["CDKN3"]
        # print DF_16
        # newDF = newDF.reindex(newDF.index.drop(0))
        return newDF.T, ID_List

    def Preprocess(self):
        LC_DBList = ['IRCR_LC14_313', 'IRCR_LC14_319', 'IRCR_LC14_320', 'IRCR_LC14_394', 'IRCR_LC14_423', 'IRCR_LC14_435', 'IRCR_LC14_436', 'IRCR_LC14_440', 'IRCR_LC14_443']
        t = stw.Timer()

        LC_train = self.Initialise_Dataset_by_Group("SGI20170607") + LC_DBList
        GBM_train = self.Initialise_Dataset_by_Group("SGI20161017") + self.Initialise_Dataset_by_Group("SGI20160613") + self.Initialise_Dataset_by_Group("MACROGEN_20151218")
        testList = self.Initialise_Dataset_by_Group("MACROGEN_20160307") + self.Initialise_Dataset_by_Group("SGI20160831") + self.Initialise_Dataset_by_Group("MACROGEN_20160516") + self.Initialise_Dataset_by_Group("TERAGEN20151119")
        labelList = self.Labeling_for_Supervised_Learning(LC_train, 0) + self.Labeling_for_Supervised_Learning(GBM_train, 1) + self.Labeling_for_Supervised_Learning(testList, 1)
        DF_train = self.Merge_DataFrame(LC_train + GBM_train + testList, 'ircr1')
        X = np.asarray(DF_train, dtype=float)
        y = np.asarray(labelList, dtype=int)
        t.stop()
        print "For Preprocessing, Time elapsed ", t.elapsed
        self.Divide_Training_Test(X, y)

    def Apply_to_non_learning_data(self):
        dataSet = self.Initialise_Dataset_by_Group("SGI20161017") + self.Initialise_Dataset_by_Group("SGI20160613") + self.Initialise_Dataset_by_Group("MACROGEN_20151218") + self.Initialise_Dataset_by_Group("MACROGEN_20160307") + \
                  self.Initialise_Dataset_by_Group("SGI20160831") + self.Initialise_Dataset_by_Group("MACROGEN_20160516") + self.Initialise_Dataset_by_Group("TERAGEN20151119")
        DF, idList = self.Merge_DataFrame(dataSet, 'ircr1')
        print "The Validation Set is  ", idList
        print "n = ", len(idList)
        self.X_val = np.asarray(DF, dtype=float)
        return idList


    def TCGA_Set(self):
        t = stw.Timer()
        with open("/home/jsgene/JK1/NGS/exec_history/TCGA_rpkm.pkl", 'r') as pp:
            DF_TCGA = pickle.load(pp)
        with open("/home/jsgene/JK1/NGS/exec_history/TCGA_array_subtype.pkl", 'r') as p:
            ID_List = pickle.load(p)
        DF_TCGA.index = ID_List
