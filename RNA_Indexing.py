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

    def RPKM_data(self, paths):
        DF = pd.DataFrame.from_csv(paths, sep='\t', header=None)
        DF = DF.reset_index()
        DF = pd.DataFrame(DF)
        SR = pd.Series(list(DF[2]), index=list(DF[0]))
        return SR

    def Merge_DataFrame(self, fileList):
        DF = []
        for afile in fileList:
            se = self.RPKM_data(afile)
            DF.append(se)
        newDF = pd.concat(DF, axis=1)
        newDF.columns = newDF.iloc[0]
        newDF = newDF.drop(["geneName"])
        # DF_16 = newDF.iloc["CDKN3"]
        # print DF_16
        # newDF = newDF.reindex(newDF.index.drop(0))
        return newDF.T

    def Preprocess(self):
        LC_train = self.Initialise_Dataset_by_Group("SGI20170607")
        GBM_train = self.Initialise_Dataset_by_Group("SGI20161017") + self.Initialise_Dataset_by_Group("SGI20160613") + self.Initialise_Dataset_by_Group("MACROGEN_20151218")
        testList = self.Initialise_Dataset_by_Group("MACROGEN_20160307") + self.Initialise_Dataset_by_Group("SGI20160831") + self.Initialise_Dataset_by_Group("MACROGEN_20160516") + self.Initialise_Dataset_by_Group("TERAGEN20151119")
        labelList = self.Labeling_for_Supervised_Learning(LC_train, 0) + self.Labeling_for_Supervised_Learning(GBM_train, 1) + self.Labeling_for_Supervised_Learning(testList, 1)
        DF_train = self.Merge_DataFrame(LC_train + GBM_train + testList )
        #DF_test = self.Merge_DataFrame(testList)
        #self.X_train = np.asarray(DF_train, dtype=float)
        #self.X_test = np.asarray(DF_test, dtype=float)
        #self.y_train = np.asarray(labelList, dtype=int)
        X = np.asarray(DF_train, dtype=float)
        y = np.asarray(labelList, dtype=int)
        from sklearn.cross_validation import train_test_split
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, random_state=123)

    def DF_to_HTML(self, newDF):
        outPath = '/data1/home/jsgene/output.html'
        DF_str = """%s""" % (newDF.to_html())
        with open(outPath, 'w') as HTML:
            HTML.write(DF_str)

##########################################################################
class visualize_RNA_expression(Setting_for_datamining, object):
    def draw_heatmap(self):
        return 0

    def draw_scatterPlot(self, X_transformed, y):
        # Plot also the training points
        colors = ['turquoise', 'darkorange']
        plt.figure(figsize=(8, 8))
        for color, i, target_name in zip(colors, [0, 1], ["LC", "GBM"]):
            plt.scatter(X_transformed[y == i, 0], X_transformed[y == i, 1], color=color, lw=2, label=target_name)
        plt.title("PCA")
        plt.legend(loc="best", shadow=False, scatterpoints=1)
        plt.axis('tight')
        plt.savefig('/home/jsgene/PCA_test.png')

    def plot_embedding(self, X, y, d, title=None):
        """Plot an embedding X with the class label y colored by the domain d."""
        x_min, x_max = np.min(X, 0), np.max(X, 0)
        X = (X - x_min) / (x_max - x_min)
        # PLOT COLORED NUMBERS
        plt.figure(figsize=(10, 10))
        ax = plt.subplot(111)
        for i in range(X.shape[0]):
            plt.text(X[i, 0], X[i, 1], str(y[i]), color=plt.cm.bwr(d[i] / 1.), fontdict={'weight': 'bold', 'size': 9})
        plt.xticks([]), plt.yticks([])  # REMOVE TICKS
        if title is not None:
            plt.title(title)
        plt.savefig('/home/jsgene/TSNE_test.png')
        
    def draw_3D_scatterPlot(self, X, y):
        from mpl_toolkits.mplot3d import Axes3D
        centres = [[1, 1], [-1, -1], [1, -1]]
        fig = plt.figure(1, figsize=(4, 3))
        plt.clf()
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
        plt.cla()

        for name, label in [('Lung_Cancer', 0), ('GBM', 1)]:
            ax.text3D(X[y == label, 0].mean(), X[y == label, 1].mean(), 200, name, horizontalalignment='center', bbox=dict(alpha=.5, edgecolor='w', facecolor='w'))
        # Reorder the labels to have colors matching the cluster results
        y = np.choose(y, [1, 2, 0]).astype(np.float)
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap=plt.cm.spectral)
        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])
        plt.savefig('/home/jsgene/PCA_test.png')

##################################################################################################

class cohort_Dataset(Setting_for_datamining, object):
    def get_GEO_SOFT(self):
        gse4271 = GEO.get_GEO(filepath="/data1/home/jsgene/KEGG/GSE4271_family.soft.gz")
        #gse1993 = GEOparse.get_GEO(filepath="/data1/home/jsgene/KEGG/GSE1993_family.soft.gz")
        #gse4422 = GEOparse.get_GEO(filepath="/data1/home/jsgene/KEGG/GSE4422_family.soft.gz")
        pivot_samples_4271 = gse4271.pivot_samples('VALUE')[3]
        print pivot_samples_4271

    def get_data_from_SQL(self, table_name, sampleID):
        import mymysql as msql
        htmlname = "/home/jsgene/sample.html"
        (con, cursor) = msql.connectDB(user='cancer', passwd='cancer', db='ircr1', host='localhost')
        DF = pd.read_sql('SELECT * FROM %s WHERE samp_id="%s"' %(table_name, sampleID), con=con)
        DF2 = pd.read_sql('SELECT DISTINCT samp_id FROM %s' % (table_name), con=con)
        print DF2
        DF_str = """%s""" % (DF2.to_html())
        with open(htmlname, 'w') as HTML:
            HTML.write(DF_str)

    def Open_RPKM_data(self, ID):
        path = self.groupName + '/%s/%s.rpkm' % (ID,ID)
        RNA_info = []
        with open(path, 'r') as f:
            lineNum = 1
            for line in f:
                (GEN, RAW_COUNTS, RPKM, ALL_READS, EXON_LENGTH) = line.split('\t')
                if lineNum == 1:
                    INDEX = [GEN, RAW_COUNTS, RPKM, ALL_READS, EXON_LENGTH]
                else:
                    Sample = [GEN, RAW_COUNTS, RPKM, ALL_READS, EXON_LENGTH]
                    Sample_info = pd.Series(data=Sample, index=INDEX)
                    RNA_info.append(Sample_info)
        DF = pd.DataFrame(RNA_info)
        return DF

###############################################################################################
class Supervised_Learning(Setting_for_datamining, object):
    def __init__(self):
        super(Supervised_Learning, self).__init__()
    def Nearest_centroid(self):
        for shrinkage in [0, .2]:
            # we create an instance of Neighbours Classifier and fit the data.
            clf = NearestCentroid(shrink_threshold=shrinkage)
            clf.fit(self.X_train, self.y_train)
            y_pred = clf.predict(self.X_test)
            test_label = []
            for i in range(len(y_pred)):
                if y_pred[i] == 0:
                    test_label.append("Lung Cancer")
                elif y_pred[i] == 1:
                    test_label.append("GBM")
            print test_label
            print(shrinkage, np.mean(self.y_train == y_pred))

        from sklearn import metrics
        accuracy = metrics.accuracy_score(self.y_test, y_pred)
        cm = metrics.confusion_matrix(self.y_test, y_pred)
        print(accuracy)
        print(cm)
        
    def k_Nearest_Neighbor(self):
        # make an instance of a k-NN classifier object
        from sklearn.neighbors import KNeighborsClassifier
        knn = KNeighborsClassifier(n_neighbors=1)
        y_pred = knn.predict(X_test)
        # calculate classification accuracy
        from sklearn import metrics
        accuracy = metrics.accuracy_score(y_test, y_pred)
        cm = metrics.confusion_matrix(y_test, y_pred)
        print(accuracy)
        print(cm)

################################################################################################

class Unsupervised_Learning(Setting_for_datamining, object):

    def __init__(self):
        super(Unsupervised_Learning, self).__init__()

    def PCA_Biopython(self):
        from Bio.Cluster import pca
        (column_mean, coordinates, Components, Eigenvalues) = pca(self.X_train)

        return Components

    def PCA_Scikit_learn(self):
        from sklearn import decomposition
        PCA = decomposition.PCA(n_components=2)
        PCA.fit(self.X_train)
        X_projection = PCA.transform(self.X_train)
        return X_projection

    def T_SNE(self):
        from sklearn.manifold import TSNE
        # RUN T-SNE
        num_test = 500
        comb_domain = np.vstack([np.tile([1., 0.], [num_test, 1]), np.tile([0., 1.], [num_test, 1])])
        tsne = TSNE(perplexity=30, n_components=2, init='pca', n_iter=3000)
        tsne_results = tsne.fit_transform(self.X_train)
        print ("T-SNE DONE.")
        plot_embedding(tsne_results, self.y_train, comb_domain.argmax(1), 'DOMAIN ADAPTATION (BLUE: SOURCE, RED: TARGET)')

if __name__ == "__main__":
    D = cohort_Dataset()
    #D.get_GEO_SOFT()
    D.get_data_from_SQL('rpkm_gene_expr', 'IRCR_GBM15_693_T03')
    #NC = Supervised_Learning()
    #NC.Preprocess()
    #NC.Nearest_centroid()
 
    #PCA = Unsupervised_Learning()
    #PCA.Preprocess()
 
    #Output = PCA.PCA_Scikit_learn()
    #PLOTS = visualize_RNA_expression()
    #PLOTS.draw_scatterPlot(Output, PCA.y_train)
