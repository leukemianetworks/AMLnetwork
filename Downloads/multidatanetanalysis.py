# multidatanetanalysis.py Version 2.0
# All mentions to GRN (Gene Regulatory Network) are referenced as
# TFG in the paper on arxiv. http://arxiv.org/ftp/arxiv/papers/1407/1407.6959.pdf
# Import Libraries
import os
import shutil
import timeit
import gzip
from zipfile import ZipFile
import numpy as np
from scipy.stats.stats import pearsonr
from scipy.stats import poisson
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Tkinter import Tk
from tkFileDialog import askopenfilename
import zipfile
import ntpath

# Define global variables
TRUE_INTERACTION = []
RESTRICT_GENE = {'tf':[],'tar':[]} 
FULL_TF = []

# Read global variable: TRUE_INTERACTION
tmp = open('TI_2486_03022014.txt').read()
for line in tmp.splitlines():
    TRUE_INTERACTION.append(line.split('\t'))
line = None
tmp = None

# Read global variable: RESTRICT_GENE
tmp = open('RG_03042014.txt').read()
for line in tmp.splitlines()[1:]:
    row = line.split('\t')
    if 'TF' in row[1]:
        RESTRICT_GENE['tf'].append(row[0])
    if 'Tar' in row[1]:
        RESTRICT_GENE['tar'].append(row[0])
row = None
line = None
tmp = None

# Read global variable: FULL_TF
tmp = open('TF_03042014.txt').read()
for line in tmp.splitlines()[1:]:
    FULL_TF.append(line.split('\t')[0])
line = None
tmp = None

# Read global variable: PPI_INTERACTION
hip = []
text = open('HIPPIE_PPI_04012014.txt').read()
for line in text.splitlines():
    tmp = line.split('\t')
    hip.append([tmp[2],tmp[5]])
HIP_INTERACTIONS = hip

# Define Error Classification
class ExpressionError(Exception):
    pass

class NetworkError(Exception):
    pass

class PipelineError(Exception):
    pass

class TIError(Exception):
    pass

def list2dict(network):
    '''
    Converts a network represented as a list to a network represented as a dictionary.
	
	Parameters:
	-----------
	network : A list
		The list-network that will be returned as a dictionary
		
	Returns:
	--------
	dict: A dictionary
		The original network in dictionary format.
		
	See Also:
	---------
	dict2list: Method to convert a network in dictionary format to list format
	
	Notes:
	------
	None.
	
	Examples:
	---------
	>>> network = readnw("path/to/a/saved/network")
	>>> networkAsDict = list2dict(network)
	
	>>> dataObject = DatasetNetworkAssemblyTool("SampleData.txt")
	>>> network = dataObject.genrnw(.1)
	>>> networkAsDict = list2dict(network)	
	
	>>> network = readnw("path/to/a/saved/network")
	>>> partialNetwork = extractnw(network, "edge", 2, 100)
	>>> networkAsDict = list2dict(partialNetwork)
    '''
	
    dict = {}
    for edge in network:
        if edge[0] not in dict.keys():
            dict[edge[0]] = [edge[1]]
        else:
            dict[edge[0]].append(edge[1])
    return dict

def dict2list(network):	
    '''
    Converts a network represented as a dictionary to a network represented as a list.
	
	Parameters:
	-----------
	network : A dictionary
		The dictionary-network that will be returned as a list.
		
	Returns:
	--------
	nw: A list
		The original network in list format.
		
	See Also:
	---------
	list2dict: Method to convert a network in dictionary format to list format
	
	Notes:
	------
	None.
	
	Examples:
	---------
	>>> network = readnw("path/to/a/saved/network")
	>>> network2 = readnw("path/to/a/saved/network2")
	>>> overlapNetwork = overlapnw(network, network2)
	>>> overlapAsList = dict2list(overlapNetwork)
    '''
    nw = []
    for key in list(network.keys()):
        if len(network[key]) > 0:
            for value in network[key]:
                nw.append([key,value])
    return nw

def readnw(nwname = "", sep = '\t'):
    '''
    Reads and returns a network file.
	
	Parameters:
	-----------
	nwname : A string
	The path to your saved network.
	The network must be in .txt format and can be zipped with a .gz extension
	
	sep : A string
	Indicates what each number will be separated by in the network file. DEFAULT
	is set to '\t' (tab).
	
	Returns:
	--------
	nw: A list
		The specified network in list format.
		
	See Also:
	---------
	extractnw: Method to get a sub-network
	
	Notes:
	------
	Network file must contain a minimum of three columns including Transcription Factor,
	Target, and method score.
	
	Examples:
	---------
	>>> network = readnw("path/to/a/saved/network.txt")
	
	>>> network2 = readnw("path/to/a/saved/network.txt.gz")
    '''
    if(nwname == ""):
	    print "Please choose the desired network to be read."
	    Tk().withdraw() 
	    nwname = askopenfilename()
    if os.path.isfile(nwname):
        if '.txt' in nwname and '.gz' in nwname:
            nwf = gzip.open(nwname,'rb')
        elif '.txt' in nwname:
            nwf = open(nwname)
        else:
            raise NetworkError('Incorrect network file format')
    else:
        raise NetworkError('Unable to find network')
		
    nw = []
    text = nwf.read()
    for line in text.splitlines():
        tmp = line.split(sep)	
        if len(tmp) < 3:
            raise NetworkError('Incorrect network data format')
        nw.append(tmp)
    return nw

def extractnw(network, ctype, start, end):
    '''
	Extract a sub-network from a given network.
	
	Parameters:
	-----------
	network : A list
	A network in list format
	
	ctype : A string
	Either 'edge' or 'percentage'
	'edge' is used when the sub-network can be obtained by indexing the passed
	in network as [start:end]
	'percentage' extracts a new network by taking a percent of the current network
    passed in -- [10:50] would take all interactions between the 10th and 50th percentile
	
	start : A string
	Indicates the starting position of the desired sub-network
	
	end : A string
	Indicates the ending position of the desired sub-network
	
	Returns:
	--------
	network: A list
		The specified sub-network in list format.
		
	Notes:
	------
	None.
	
	Examples:
	---------
	>>> network = readnw("path/to/a/saved/network.txt")
	>>> newNetwork = extractnw(network, "percentage", 1, 100)
	
	>>> network2 = readnw("path/to/a/saved/network.txt.gz")
	>>> newNetwork2 = extractnw(network2, "edge", 2, 435)
    '''
    if ctype == 'edge':
        start = start - 1
        return network[start:end]
    if ctype == 'percentage':
        return network[int(round(start/100.0*len(network))):int(round(end/100.0*len(network)))]
	
def overlapnw(nw1, nw2):
    '''
    Returns the dictionary containing the overlaps between two input networks
	and also the count of such overlapping.
	
	Parameters:
	-----------
	nw1 : A list
	A network in list format
	
	nw2 : A list
	A network in list format

	Returns:
	--------
	overlap: A dictionary
		The key of the dictionary is each
		
	count: An integer
		The total number of overlaps	

	
	Examples:
	---------
	>>> network1 = readnw("path/to/a/saved/network1.txt")
	>>> network2 = readnw("path/to/a/saved/network2.txt")
	>>> overlapNetwork = overlapnw(network1, network2)
    '''
    if type(nw1) == list:
        nw1 = list2dict(nw1)
    if type(nw2) == list:
        nw2 = list2dict(nw2)
    overlap = {}
    count = 0
    for key in list(nw1.keys()):            
        if key in nw2.keys():
            ol = set(nw1[key]).intersection(set(nw2[key]))
            overlap[key] = ol
            count += len(ol)
    return overlap, count

def cvcut(gene, expr, cutoff):
    '''
    Filtering expression with a given CV. Called by optcv in the
	DatasetNetworkAssemblyTool Class
	
	Parameters:
	-----------
	gene : A list
	A list of genes
	
	expr : A list
	A network in list format
	
	cutoff : A float
	The desired cutoff number

	Returns:
	--------
	genecut: A list
		The new list of genes after the cutoff is implemented
		
	exprcut: An integer
		The expr list updated from the cut
	Notes:
	------
	It is calculated by the standard deviation divided by the mean.
	
	Examples:
	---------
	>>> gene, expr = cvcut(self.rgene, self.rexpr, self.cv)
    '''
    genef = []
    exprf = []
    for g in range(len(gene)):
        cv = np.std(expr[:,g])/np.mean(expr[:,g])
        if cv > cutoff:
            genef.append(gene[g])
            exprf.append(expr[:,g])
    return genef, np.transpose(exprf)

def runpcorr(tfs, gene, expr, verbal = True):
    '''
    Generate network via Pearson Correlation
	
	Parameters:
	-----------
	tfs : A list
	A list of transcription factors
	
	gene : A list
	A list of genes
	
	expr : A float
	The desired cutoff number
	
	verbal : A boolean
	Indicates whether program will show updates on how far the algorithm has
	progressed through the list of genes

	Returns:
	--------
	nw : A list
	The network generated
	
	Notes:
	------
	P corr relates to Pearson correlation.
	
	Examples:
	---------
	Primarily called when generating networks inside of the DatasetNetworkAssemblingTool
	class. Specifically, runpcorr() is called in optcv(), genfnw(), and genrnw().

    '''
    # Network Generation
    nw = []
    count = 1
    for tar in gene:
        if verbal:
            print 'Processing Gene %i/%i' % (count,len(gene))
        start = timeit.default_timer()
        for tf in tfs:
            if tf == tar or tf not in gene:
                continue
            else:
                coeff, p = pearsonr(expr[:,gene.index(tf)],expr[:,gene.index(tar)])
                nw.append([tf,tar,coeff])
        stop = timeit.default_timer()
        if verbal:
            print '\tProcess time: %s' % str(stop-start)
        count += 1
    nw.sort(key=lambda ind:ind[2], reverse=True)
    return nw

class SingleDatasetNetworkAssembling:
    '''
    SingleDatasetNetworkAssembling reads in Microarray or RNASeq expression data,
    and generates a GRN and PPI subnetwork. The dataset is first restricted
    into a smaller size data. A set of benchmarked
	"True Interactions"(TI)(Experimentally/literature verified interaction
	is used for both coefficient of variation (cv) filtering optimization and
    statistical analysis. CV is optimized base on the overlap between top 100
	ranked interactions of the predicted GRN subnetwork from restricted expression
	profile and TI. Statistical analysis is also performed using the optimized CV.
    With the optimized CV, a full size GRN subnetwork is be generated Using poission
    distribution, a statistically significant correlation of coefficient (cc) is reported.
    The filtered GRN and PPI subnetwork is inferred base on the cc.
    '''
    
    def __init__(self, expfile = "", sep = '\t', form = 'row.as.condition'):
        '''
		Constructor for SingleDatasetNetworkAssembling.
        
		Parameters:
		-----------
        expfile: A file
		Expression data is suggested to be processed. Suggested
		processing as follow:
            For Microarray data:
                Data transformation: Log2 Transformation
                Background adjustment: RMA2
                Normalization: Quantile
                Summarization: Median Polish
                Multi-probe to single-gene mapping: Median
            For RNASeq data:
                Normalization: RPKM
                Data transformation: Ln(RPKM+1)
				
        sep: A String, optional
		Default separator is tab-delimited file
		
		form: A boolean, optional
		Expression data format can be either
		'row.as.condition' or 'row.as.gene'. Default is 'row.as.condition', which
		is the format used in the following procedures.
        
		Returns:
		--------
		object: SingleDatasetNetworkAssembling 
		An instance of the class
		
		Notes:
		------
		None.
		
		Examples:
		---------
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/data.txt", compr = False)
		
		>>> classInstance2 = SingleDatasetNetworkAssembling("path/to/compressed/data.txt.gz")
		
		>>> classInstance3 = SingleDatasetNetworkAssembling("path/to/row.as.gene/data.txt", compr = False, form = "row.as.gene")
        '''
        
        self.expr = []
        self.gene = []
        self.cond = []
        self.rexpr = []
        self.rgene = []
        self.plname = ''
        
        self.cv = None
        self.cc = None
        self.rgrn = None
        self.fgrn = None
        self.grn = None
        self.fppi = None
        self.ppi = None
        self.expf = None
		
        if(expfile == ""):
            print "Please choose the desired file to be read."
            Tk().withdraw() 
            expfile = askopenfilename()
        
        # Check compression and file existence
        # Read expression
        print(expfile)
        expfile = ntpath.basename(expfile)
        if '.gz' in expfile:
			if os.path.isfile(expfile):
				exprf = gzip.open(expfile,'rb')
				indexOfFileExtension = expfile.rfind(".", 0, expfile.rfind("."))
				plname = expfile[0:indexOfFileExtension]
			else:
				raise ExpressionError('Unable to find expression data.')
        elif '.zip' in expfile:
            fh = open(expfile, "rb")
            zipObject = ZipFile(fh)
            for name in zipObject.namelist():
                outpath = os.path.dirname(os.path.realpath("multidatanetanalysis.py"))
                zipObject.extract(name, outpath)
            fh.close()
            print("Please choose the correct file from the extracted files")
            Tk().withdraw() 
            expfile = askopenfilename()
            indexOfFileExtension = expfile.rfind(".")
            plname = expfile[0:indexOfFileExtension]
            
        if os.path.isfile(expfile):
            exprf = open(expfile,'r')
            indexOfFileExtension = expfile.rfind(".")
            plname = expfile[0:indexOfFileExtension]
            if expfile[indexOfFileExtension:] not in (".txt", ".gz", ".zip"):
                raise ExpressionError("Only support .gz, .zip, and .txt file extensions")
            
        else:
            raise ExpressionError('Unable to find expression data.')
        
        compr = True
        if expfile[indexOfFileExtension + 1:] == ".txt":
            print(expfile[indexOfFileExtension:])
            compr = False
        # Read data
        data = exprf.read()
        expr = []
        gene = []
        cond = []
        if form == 'row.as.condition':
            gene = data.splitlines()[0].split(sep)[1:]
        else:
            cond = data.splitlines()[0].split(sep)[1:]
        for line in data.splitlines()[1:]:
            tmp = line.split(sep)
            if form == 'row.as.condition':
                cond.append(tmp[0])
            else:
                gene.append(tmp[0])
            expr.append(tmp[1:])
        expr = np.array(expr).astype(float)
        data = None
        tmp = None
        if form != 'row.as.condition':
            expr = np.transpose(expr)
        
        self.expr = expr
        self.gene = gene
        self.cond = cond
        
        # Check data
        if np.size(expr,0) != len(cond) or np.size(expr,1) != len(gene):
            raise ExpressionError('Incorrect expression format.')
        
        # Restrict data
        tigene = list(set(RESTRICT_GENE['tf'] + RESTRICT_GENE['tar']))
        for i in range(len(gene)):
            if gene[i] in tigene:
                self.rgene.append(gene[i])
                self.rexpr.append(expr[:,i])
        self.rexpr = np.transpose(np.array(self.rexpr).astype(float))
        
        # Create pipeline folder in current directory
        self.plname = '%s' % plname
        if os.path.isdir(self.plname):
            (self.cv,self.cc,self.rgrn,self.fgrn,self.grn,self.fppi,self.ppi, self.expf) = self.getstatus()
        else:
            os.mkdir(self.plname)
            # Create and keep track of pipeline status
            self._UpdateStatus()
        
        # Check if copy expression data to pipeline folder
        if self.expf == None:
            shutil.copy(expfile,'%s/%s' % (self.plname,expfile))		
            if compr:
                self.expf = '%s/%s' % (self.plname,expfile)
            else:
                tmp = open('%s/%s' % (self.plname,expfile)).read()
                gzip.open('%s/%s.gz' % (self.plname,expfile), 'wb').write(tmp)
                os.remove('%s/%s' % (self.plname,expfile))
                self.expf = '%s/%s.gz' % (self.plname,expfile)
    
    def _UpdateStatus(self):
        '''
        Updates the SingleDatasetNetworkAssembling object status.
        
        Parameters:
        -----------
        None.

        Returns:
        --------
        None:
        
        Notes:
        ------   
        The generated status file will be output into a file named current/directory/YYY/Status.txt
        where YYY indicates the name of the data set being generating.
        
        Examples:
        ---------
        >>> classInstance = SingleDatasetNetworkAssembling("path/to/data/txt.gz")
        >>> classInstance.getstatus()
        '''
        
        log = ''
        log += 'Pipeline: %s' % self.plname
        log += '\n\nExpression profile: %s' % self.plname.replace('PL_','')
        log += '\nExpression profile location: %s' % self.expf
        log += '\nExpression properties:'
        if self.cv == None:
            log += '\n\t-# of genes: %s' % len(self.gene)
        else:
            gene, expr = cvcut(self.rgene, self.rexpr, self.cv)
            log += '\n\t-# of genes: %s (%i genes remain)' % (len(self.gene),len(gene))
        log += '\n\t-# of conditions: %s' % len(self.cond)
        log += '\n\nRestricted Expression properties:'
        log += '\n\t-# of genes: %s' % len(self.rgene)
        log += '\n\t-# of conditions: %s' % len(self.cond)
        log += '\n\nPipeline network status:'
        log += '\n\t-CV Optimization: %s' % str(self.cv)
        log += '\n\t-Significant CC value: %s' % str(self.cc)
        log += '\n\t-CV Opt Restricted GRN subnetwork: %s' % str(self.rgrn)
        log += '\n\t-CV Opt Full GRN subnetwork: %s' % str(self.fgrn)
        log += '\n\t-CV Opt CC Cut complete GRN subnetwork: %s' % str(self.grn)
        log += '\n\t-Full PPI subnetwork: %s' % str(self.fppi)
        log += '\n\t-CC Cut complete PPI subnetwork: %s' % str(self.ppi)
        open('%s/Status.txt' % self.plname, 'w').write(log)
    
    def status(self):
        '''
        print the SingleDatasetNetworkAssembling object history.
        
        Parameters:
        -----------
        None.

        Returns:
        --------
        None:
        
        Notes:
        ------   
        None
        
        Examples:
        ---------
        >>> classInstance = SingleDatasetNetworkAssembling("path/to/data/txt.gz")
        >>> classInstance.getstatus()
        '''
        
        self._UpdateStatus()
        print open('%s/Status.txt' % self.plname).read()
    
    def getstatus(self):
        tmp = open('%s/Status.txt' % self.plname).read().splitlines()
        expf = tmp[2].split(': ')[1]
        if not os.path.isfile(expf):
            expf = None
        cv = eval(tmp[13].split(': ')[1])
        cc = eval(tmp[14].split(': ')[1])
        rgrn = tmp[15].split(': ')[1]
        if not os.path.isfile(rgrn):
            rgrn = None
        fgrn = tmp[16].split(': ')[1]
        if not os.path.isfile(fgrn):
            fgrn = None
        grn = tmp[17].split(': ')[1]
        if not os.path.isfile(grn):
            grn = None
        fppi = tmp[18].split(': ')[1]
        if not os.path.isfile(fppi):
            fppi = None
        ppi = tmp[19].split(': ')[1]
        if not os.path.isfile(ppi):
            ppi = None
        return (cv, cc, rgrn, fgrn, grn, fppi, ppi, expf)
    
    def optcv(self):
        '''
        Updates the SingleDatasetNetworkAssembling object as well as creates
        a plot based on the number of overlaps in the network and their respective 
        CV value.
        
        Parameters:
        -----------
        None.

        Returns:
        --------
        None:
        
        Notes:
        ------   
        CV optimization will be done in a range of CV = [0,1.0] based on the
        overlap between top 100 ranked interactions
        of the predicted GRN subnetwork from restricted expression profile and TI. The
        distribution of overlap against CV is expected to be normally distributed.
        Hence, there is a maximum CV which optimize the overlapping. If there are
        multiple CV which optimize the overlapping, median CV of the lowest interval will be used.
        For more information, reference
        http://arxiv.org/ftp/arxiv/papers/1407/1407.6959.pdf
        
        Examples:
        ---------
        >>> classInstance = SingleDatasetNetworkAssembling("path/to/data/txt.gz")
        >>> classInstance.optcv()
        '''
        
        print 'Initiating CV Optimization...'
        tfs = RESTRICT_GENE['tf']
        pp = PdfPages('%s/CV_Optimization.pdf' % self.plname)
        text = ''
        
        # First iteration
        print '1st iteration [0.0,1.0] with 0.1 interval'
        cvs = np.arange(0.0,1.1,0.1)
        ols = []
        text += 'First iteration\t' + '\t'.join(list(cvs.astype(str))) + '\nOverlap'
        for cv in cvs:
            print 'Processing %s CV' % str(cv)
            gene, expr = cvcut(self.rgene, self.rexpr, cv)
            nw = runpcorr(tfs, gene, expr, verbal=False)
            nw = [[i,j] for i,j,k in nw]
            nw = extractnw(nw, 'edge', 1, 100)
            if len(nw) < 100:
                ols.append(0)
                text += '\t0'
            else:
                ol, count = overlapnw(nw, TRUE_INTERACTION)
                ols.append(count)
                text += '\t%i' % count
        mind = ols.index(np.max(ols))
        tmp = [mind]
        for ind in range(mind,len(cvs)):
            if ols[ind] == ols[mind]:
                tmp.append(ind)
        mind = int(round(np.median(tmp)))
        plt.scatter(cvs[:mind], ols[:mind], marker='o', c='k')
        plt.scatter(cvs[mind], ols[mind], marker='o', c='r')
        plt.scatter(cvs[mind+1:], ols[mind+1:], marker='o', c='k')
        plt.xlabel('CVs')
        plt.ylabel('# of overlapping interactions')
        plt.title('1st iteration CV Optimization')
        plt.savefig(pp,format='pdf')
        plt.close()
        
        # Second iteration
        if mind == 0:
            print '2nd iteration [0.00,%.2f] with 0.01 interval' % (cvs[mind+1])
            cvs = np.arange(0, cvs[mind+1]+0.01,0.01)
        else:
            print '2nd iteration [%.2f,%.2f] with 0.01 interval' % (cvs[mind-1],cvs[mind+1])
            cvs = np.arange(cvs[mind-1], cvs[mind+1]+0.01,0.01)
        ols = []
        text += '\nSecond iteration\t' + '\t'.join(list(cvs.astype(str))) + '\nOverlap'
        for cv in cvs:
            print 'Processing %s CV' % str(cv)
            gene, expr = cvcut(self.rgene, self.rexpr, cv)
            nw = runpcorr(tfs, gene, expr, verbal=False)
            nw = [[i,j] for i,j,k in nw]
            nw = extractnw(nw, 'edge', 1, 100)
            if len(nw) < 100:
                ols.append(0)
                text += '\t0'
            else:
                ol, count = overlapnw(nw, TRUE_INTERACTION)
                ols.append(count)
                text += '\t%i' % count
        mind = ols.index(np.max(ols))
        tmp = [mind]
        for ind in range(mind,len(cvs)):
            if ols[ind] == ols[mind]:
                tmp.append(ind)
        mind = int(round(np.median(tmp)))
        plt.scatter(cvs[:mind], ols[:mind], marker='o', c='k')
        plt.scatter(cvs[mind], ols[mind], marker='o', c='r')
        plt.scatter(cvs[mind+1:], ols[mind+1:], marker='o', c='k')
        plt.xlabel('CVs')
        plt.ylabel('# of overlapping interactions')
        plt.title('2nd iteration CV Optimization')
        plt.savefig(pp,format='pdf')
        plt.close()
        
        # Optimized CV
        open('%s/CV_Optimization.txt' % self.plname,'w').write(text)
        pp.close()
        self.cv = cvs[mind]
        
        # Update pipeline history
        self._UpdateStatus()
    
    def genrestrictgrn(self, cv, output = False, compr = True):
        '''
        Generates and returns a restricted version of GRN subnetwork based on interactions 
		in the restricted expression profile via Pearson Correlation.
	
		Parameters:
		-----------
		cv : A float
		Represents the coefficient of variance
		
		output : A boolean, optional
		Represents if you would like the output saved to a file.
		Filename will be in in the folder %s where %s represents your data
		set name. The file itself will be named Restricted_GRN_Subnetwork_XXcv where
		'XX' represents the cv the network was generated with.
		Default is set to False.
		
		compr : A boolean, optional
		Indicates whether the outputted file will be compressed. Default is set
		to True.
		
		Returns:
		--------
		rgrn : A list
		The restricted GRN subnetwork generated
		
		Notes:
		------
		None.
		
		Examples:
		---------
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.genrestrictgrn(.09)
		
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.genrestrictgrn(.09, output = True, compr = False)
        '''
        
        rgrn = []
        tfs = RESTRICT_GENE['tf']
        
        print 'Generating restricted GRN interactions via Pearson Correlation...'
        gene, expr = cvcut(self.rgene, self.rexpr, cv)
        rgrn = runpcorr(tfs, gene, expr, verbal=False)
        
        if output:
            result = []
            for row in np.array(rgrn).astype(str):
                result.append('\t'.join(row))
            if compr:
                gzip.open('%s/Restricted_GRN_Subnetwork_%.2fcv.txt.gz' % (self.plname, cv),'wb').write('\n'.join(result))
            else:
                open('%s/Restricted_GRN_Subnetwork_%.2fcv.txt' % (self.plname, cv),'w').write('\n'.join(result))
        return rgrn
    
    def anlzrestrictgrn(self, cv = 'optcv', binsize = 100, verbal = True):
        '''
        Analyses restricted GRN subnetwork against True Interactions and creates a
		TI Cumlative Overlapping and TI Binned Overlapping Distribution plots
	
		Parameters:
		-----------
		cv : A string/float, optional
		Default is set to 'optcv' which gets data's optimal cv via either self.cv or optcv()
		if it has yet too be initialized. User can also input other cv to be used.
		
		binsize : An integer, optional
		Default is set to 100.
		
		verbal: A boolean, optional
		Default set to True. When True, algorithm outputs which TI overlap
		bin it is currently processing.
		
		Returns:
		--------
		None
		
		Notes:
		------
		Plots can be found in the XXX folder where 'XXX' represents your data
		file name.
		Function creates a bar graph showing the number of TI Hits in each bin.
        The size of each bin corresponds to the different intervals of the 
        correlation coefficient.
		
		Examples:
		---------
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.anlzrestrictgrn()
		
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.anlzrestrictgrn(.09, binsize = 200, verbal = False)
        '''
        
        if cv == 'optcv':
            if self.cv == None:
                self.optcv()
            rgrn = self.genrestrictgrn(self.cv, output = True)
            self.rgrn = '%s/Restricted_GRN_Subnetwork_%.2fcv.txt.gz' % (self.plname, self.cv)
        else:
            rgrn = self.genrestrictgrn(cv, output = True)
            
        bins = len(rgrn)/binsize
        nw = extractnw(rgrn, 'edge', 1, binsize*bins)
        tidic = list2dict(TRUE_INTERACTION)
        
        # Process TI Overlap and slope for each bin. Also calculate skewness of the corresponding TI Overlap distribution
        overlaps = np.array([[0] * (len(tidic.keys())+1)] * bins)
        overlapd = []
        slopes = []
        counts = []
        text = 'TF_Symbol\tTI_Target#'
        for bin in range(bins):
            p1 = float(bin)/bins * 100
            start = int(round(len(nw)*p1/100))
            p2 = (bin+1.0)/bins * 100
            end = int(round(len(nw)*p2/100))
            text += '\t%s%%(%i)_Edges' % (str(p2),end)
            if verbal:
                print 'Processing TI Overlap bin: %s%% to %s%%...' % (str(p1),str(p2))
            tmpnw = extractnw(nw, 'edge', start, end)
            overlap, count = overlapnw(tidic, tmpnw)
            overlaptf = [0] * len(tidic.keys())
            for key in list(overlap.keys()):
                keyindex = list(tidic.keys()).index(key)
                overlaptf[keyindex] =len(overlap[key])
            if bin == 0:
                overlaps[bin,] = [count] + overlaptf
            else:
                overlaps[bin,] = overlaps[bin-1,] + ([count]+overlaptf)
            slopes.append(float(count)/(end-start))
            counts.append(count)
            overlapd += [p1] * count
        overlaps = np.array(overlaps)
        overlaps = np.transpose(overlaps)
        mean = np.mean(overlapd)
        median = np.median(overlapd)
        std = np.std(overlapd)
        skew = 3 * (mean - median) / std
        text += '\nCumlative Overlap\t%i' % len(TRUE_INTERACTION)
        for i in range(bins):
            text += '\t%i' % overlaps[0,i]
        text += '\nBinned Overlap\tNA'
        for i in range(bins):
            text += '\t%i' % counts[i]
        text += '\nBinned Overlap:Interactions\tNA'
        for i in range(bins):
            text += '\t%.4f' % slopes[i]
        text += '\nBinned Overlap:TransFac\tNA'
        for i in range(bins):
            text += '\t%.4f' % (float(counts[i])/len(TRUE_INTERACTION))
        keyindex = 1
        for key in list(tidic.keys()):
            text += '\n%s\t%i' % (key,len(tidic[key]))
            for i in range(bins):
                text += '\t%i' % overlaps[keyindex,i]
            keyindex += 1
        
        # Write TI Overlap result matrix.
        open('%s/Restricted_GRN_Subnetwork_%.2fcv_TIOverlapMatrix_%iBinsize.txt' % (self.plname,self.cv,binsize),'w').write(text)
        
        # Plot binned TI Overlap distribution(bar graph), and cumlative TI Overlap(line graph)
        pp = PdfPages('%s/Restricted_GRN_Subnetwork_%.2fcv_TIOverlapGraph_%iBinsize.pdf' % (self.plname,self.cv,binsize))
        
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(np.array(range(bins+1))/(bins/100.0),[0]+list(overlaps[0,]),'k-')
        ax.set_xlabel('Cumulative % Interactions')
        ax.set_ylabel('# of Overlapping Interactions')
        ax.set_title('TI Cumlative Overlapping')
        fig.savefig(pp,format='pdf')
        plt.close()
        
        n, num_bins, patches = plt.hist(overlapd, bins, range=(0,100), facecolor='green', alpha=0.5)
        plt.xlabel('% Interactions')
        plt.ylabel('# of Overlapping Interactions')
        plt.title('TI Binned Overlapping Distribution (Skewness = %.3f)' % skew)
        # Tweak spacing to prevent clipping of ylabel
        plt.subplots_adjust(left=0.15)
        plt.savefig(pp,format='pdf')
        plt.close()
        
        pp.close()
        
        # Update pipeline history
        self._UpdateStatus()
    
    def genfullgrn(self, cv = 'optcv', compr = True):
        '''
        Generates a full size GRN subnetwork with the optimized CV.
	
		Parameters:
		-----------
		cv : A string/float, optional
		Default is set to 'optcv' which gets data's optimal cv via either self.cv or optcv()
		if self.cv has yet too be initialized. User can also input other cv to be used.
		
		compr : A boolean, optional
		Default set to True. Indicates whether the output network will be saved
		in a compressed file (gz format)
		
		Returns:
		--------
		None
		
		Notes:
		------
		The generated GRN subnetwork will be output into a file named current/directory/YYY/GRN_Subnetwork_%XXXcv
		where YYY indicates the name of the data set being generating and XXX indicates
		the cv used in generating the network.
		
		Examples:
		---------
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.genfullgrn()
		
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.genfullgrn(.09, compr = False)
        '''
        
        if cv == 'optcv':
            if self.cv == None:
                self.optcv()
            gene, expr = cvcut(self.gene, self.expr, self.cv)
        else:
            gene, expr = cvcut(self.gene, self.expr, cv)
        
        fgrn = []
        tfs = FULL_TF
        
        print 'Generating full GRN interactions via Pearson Correlation...'
        fgrn = runpcorr(tfs, gene, expr, verbal=True)
        
        result = []
        for row in np.array(fgrn).astype(str):
            result.append('\t'.join(row))
        if compr:
            gzip.open('%s/Full_GRN_Subnetwork_%.2fcv.txt.gz' % (self.plname, self.cv),'wb').write('\n'.join(result))
            if cv == 'optcv':
                self.fgrn = '%s/Full_GRN_Subnetwork_%.2fcv.txt.gz' % (self.plname, self.cv)
        else:
            open('%s/Full_GRN_Subnetwork_%.2fcv.txt' % (self.plname, cv),'w').write('\n'.join(result))
            if cv == 'optcv':
                self.fgrn = '%s/Full_GRN_Subnetwork_%.2fcv.txt' % (self.plname, self.cv)
        
        # Update pipeline history
        self._UpdateStatus()

    def genfullppi(self, cv = 'optcv', compr = True):
        '''
        Generates a full size PPI subnetwork with the optimized CV.
	
		Parameters:
		-----------
		cv : A string/float, optional
		Default is set to 'optcv' which gets data's optimal cv via either self.cv or optcv()
		if self.cv has yet too be initialized. User can also input other cv to be used.
		
		compr : A boolean, optional
		Default set to True. Indicates whether the output network will be saved
		in a compressed file (gz format)
		
		Returns:
		--------
		None
		
		Notes:
		------
		The generated PPI subnetwork will be output into a file named current/directory/YYY/PPI_Subnetwork_%XXXcv
		where YYY indicates the name of the data set being generating and XXX indicates
		the cv used in generating the network.
		
		Examples:
		---------
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.genfullppi()
		
		>>> classInstance = SingleDatasetNetworkAssembling("path/to/your/data.txt.gz")
		>>> classInstance.genfullppi(.09, compr = False)
        '''
        
        gene = self.gene
        expr = self.expr
        fppi = []
        hip = HIP_INTERACTIONS
        used = []
        
        print 'Generating full PPI interactions via Pearson Correlation...'
        cur = 0
        for i in range(len(hip)):
            if float(i)/len(hip)*100 > cur:
                print 'Processing PPI %i%% to %i%%' % (cur,cur+1)
                cur += 1
            ppi = hip[i]
            if ppi[0] in gene and ppi[1] in gene and ppi[0] != ppi[1] and ppi not in used:
                coeff, p = pearsonr(expr[:,gene.index(ppi[0])],expr[:,gene.index(ppi[1])])
                fppi.append([ppi[0],ppi[1],coeff])
                used.append(ppi)
        
        result = []
        for row in np.array(fppi).astype(str):
            result.append('\t'.join(row))
        if compr:
            gzip.open('%s/Full_PPI_Subnetwork_%.2fcv.txt.gz' % (self.plname, self.cv),'wb').write('\n'.join(result))
            if cv == 'optcv':
                self.fppi = '%s/Full_PPI_Subnetwork_%.2fcv.txt.gz' % (self.plname, self.cv)
        else:
            open('%s/Full_PPI_Subnetwork_%.2fcv.txt' % (self.plname, cv),'w').write('\n'.join(result))
            if cv == 'optcv':
                self.fppi = '%s/Full_GRN_Subnetwork_%.2fcv.txt' % (self.plname, self.cv)
        
        # Update pipeline history
        self._UpdateStatus()
    
    def poissoncut(self, fgrnfile = None, fppifile = None, pval = 0.10, rate = 2.0):
        '''
        Extract significant interactions for both GRN and PPI subnetwork base on given p-value using poisson
		distribution. This function also update the corresponding cc to the SingleDatasetNetworkAssembling object.
	    
		Parameters:
		-----------
		fgrnfile : A String, optional
		Default is set to None. When none, function generates full GRN subnetwork via 
		genfullgrn(). User also has option to input the directory of a saved network
		to run the calculation on.
        
        fppifile : A String, optional
		Default is set to None. When none, function generates full PPI subnetwork via 
		genfullppi(). User also has option to input the directory of a saved network
		to run the calculation on.
		
		pval : A float, optional
		Default set to .05. Indicates the p-value to be used during the poissoncut.
		
		rate : An float, optional
		Default set to 1.0. The rate of TI hit.
		
		Returns:
		--------
		None
		
		Notes:
		------
		Will output data in file named YYY/Network_XXXfcv_XXXfpval where YYY
		represents the name of the current data set being analyzed and YYY
		represents the cv used for the data.
		
		Examples:
		---------
		>>> classInstance = DatasetNetworkAssemblingTool("path/to/your/data.txt.gz")
		>>> classInstance.poissoncut()
		
		>>> classInstance = DatasetNetworkAssemblingTool("path/to/your/data.txt.gz")
		>>> classInstance.poissoncut("/path/to/saved/fullNetwork")
		
		>>> classInstance = DatasetNetworkAssemblingTool("path/to/your/data.txt.gz")
		>>> classInstance.poissoncut("/path/to/saved/fullNetwork", pval = .10, bins = 100, binsize = 7000)
        '''
        
        if self.cv == None:
            self.optcv()
        
        print 'Loading GRN Subnetwork...'
        if fgrnfile == None:
            if self.fgrn == None:
                self.genfullgrn()
            fgrn = readnw(self.fgrn)
        else:
            if os.path.isfile(fgrnfile):
                self.fgrn = fgrnfile
                fgrn = readnw(fgrnfile)
            else:
                raise NetworkError('Unable to find full size GRN subnetwork.')
        
        print 'Loading PPI Subnetwork...'
        if fppifile == None:
            if self.fppi == None:
                self.genfullppi()
            fppi = readnw(self.fppi)
        else:
            if os.path.isfile(fppifile):
                self.fppi = fppifile
                fppi = readnw(fppifile)
            else:
                raise NetworkError('Unable to find full size PPI subnetwork.')
        
        fgrnc = [k for i,j,k in fgrn]
        fgrn = [[i,j] for i,j,k in fgrn]
        ti = [[i,j] for i,j in TRUE_INTERACTION if [i,j] in fgrn]
        tfs = FULL_TF
        ftfs = set(tfs).intersection(set(self.gene))
        gene, expr = cvcut(self.gene, self.expr, self.cv)
        fnsize = len(fgrn)
        binsize = int(rate/(float(len(ti))/fnsize))
        bins = int(fnsize/binsize)
        
        # Calculate Poisson Distribution Cutoff
        # Also update the corresponding cc and plot TI hit histogram
        print 'Processing Poisson cutoff...'
        cutoff = poisson.ppf((1.0-pval),float(len(ti))/fnsize*binsize)
        h = []
        e = []
        grn = None
        for i in range(bins):
            tmp = extractnw(fgrn, 'edge', (i*binsize)+1, (i+1)*binsize)
            overlaps, count = overlapnw(tmp,ti)
            h.append(count)
            e.append(i)
            if grn == None and count < cutoff:
                grn = extractnw(fgrn, 'edge', 1, i*binsize)
                vl = i
            if grn != None and len(h) >= 50:
                break
        # Process Poisson filtered GRN subnetwork
        print 'Processing Poisson cutoff GRN subnetwork...'
        result = ['Rank\tTranscription Factor\tTarget\tScore(Pipeline Specific)\tif TransFac']
        for i in range(len(grn)):
            if grn[i] in ti:
                result.append('%i\t%s\t%s\t%s' % (i+1,'\t'.join(grn[i]),fgrnc[i],'1'))
            else:
                result.append('%i\t%s\t%s\t%s' % (i+1,'\t'.join(grn[i]),fgrnc[i],'0'))
        self.cc = round(float(fgrnc[i]),2)
        
        # Plot TI Hit histogram
        pp = PdfPages('%s/GRN_TIHit_%.2fcv_%.2fcc_%.2fpval_%.1frate.pdf' % (self.plname,self.cv,self.cc,pval,rate))
        r = plt.bar(e, h, color='b')
        plt.xlabel('Bins. Bin size = %i. Bins = %i' % (binsize,len(h)))
        plt.ylabel('TI Hit')
        plt.title('TI Hit of binned interactions. Avg rate = %s' % str(rate))
        plt.axvline(x=vl,color='r')
        plt.savefig(pp,format='pdf')
        plt.close()
        pp.close()
        
        # Write Poisson filtered GRN subnetwork
        open('%s/GRN_Subnetwork_%.2fcv_%.2fcc_%.2fpval_%.1frate.txt' % (self.plname,self.cv,self.cc,pval,rate),'w').write('\n'.join(result))
        self.grn = '%s/GRN_Subnetwork_%.2fcv_%.2fcc_%.2fpval_%.1frate.txt' % (self.plname,self.cv,self.cc,pval,rate)
        
        # Process Poisson filtered PPI subnetwork
        print 'Processing Poisson cutoff PPI subnetwork...'
        result = ['Rank\tGene1\tGene2\tScore(Pipeline Specific)']
        for i in range(len(fppi)):
            if float(fppi[i][2]) > self.cc:
                result.append('%i\t%s' % (i+1,'\t'.join(fppi[i])))
        
        # Write Poisson filtered PPI subnetwork
        open('%s/PPI_Subnetwork_%.2fcv_%.2fcc_%.2fpval_%.1frate.txt' % (self.plname,self.cv,self.cc,pval,rate),'w').write('\n'.join(result))
        self.ppi = '%s/PPI_Subnetwork_%.2fcv_%.2fcc_%.2fpval_%.1frate.txt' % (self.plname,self.cv,self.cc,pval,rate)
    
        # Update pipeline history
        self._UpdateStatus()


class MultiDatasetNetworkAssembling:
    ''' MultiDatasetNetworkAssembling allows for multiple network files to be compared.
        Overlap networks can be created by using the intppi() and intgrn() methods once
        all desired expression files have been analyzed.
        '''
    
    def __init__(self, expfiles = [], numFiles = -1, sep = '\t', form = 'row.as.condition'):
        '''
		Constructor for MultiDatasetNetworkAssembling.
        
		Parameters:
		-----------
        expfiles: A list, optional
            A list of locations of expfiles desired to be analyzed. If none are given
            user will be prompted to pick a file 'numFiles' times. Default set to empty
            list.
            
        numFiles: an integer, optional
            The number of files the user wishes to analyzed
            
        sep: A String, optional
		Default separator is tab-delimited file
		
		form: A boolean, optional
		Expression data format can be either
		'row.as.condition' or 'row.as.gene'. Default is 'row.as.condition', which
		is the format used in the following procedures.
        
		Returns:
		--------
		object: MultiDatasetNetworkAssembling 
		An instance of the class
		
		Notes:
		------
		If expression file is found not to be previously analyzed, function will
        automatically call SingleDatasetNetworkAssembling() on the respective file.
		
		Examples:
		---------
		>>> classInstance = MultiDatasetNetworkAssembling(numFiles = 4)
		
		>>> classInstance2 = MultiDatasetNetworkAssembling(["path/to/compressed/data.txt.gz", "path/to/compressed/data2.txt.gz"])
        
        '''
        self.exps = []
        self.pls = []
        self.grns = []
        self.ppis = []
        
        self.igrn = None
        self.ippi = None
        
        if numFiles == -1 and expfiles == []:
            raise ValueError("Either need to set expfiles manually or numFiles to the number of expFiles to be compared")
        print expfiles
        if(expfiles == []):
            for i in range(0,numFiles) :
                print "Please choose the desired file to be read."
                Tk().withdraw() 
                expfile = askopenfilename()
                if expfile == "":
                    expfiles = []
                    raise ValueError("expfile file not chosen")
                expfiles.append(expfile)
                print expfile
        
                
        exp_len = len(expfiles)
        for exp in range (exp_len):
            pl = 'PL_%s' % exp
            print 'Processing pipeline: %s' % pl
            self.pls.append(pl)
            sdna = SingleDatasetNetworkAssembling(expfiles[0])
            if sdna.grn == None or sdna.ppi == None:
                SingleDatasetNetworkAssembling.poissoncut(sdna)
            self.grns.append(sdna.grn)
            self.ppis.append(sdna.ppi)
            expfiles.pop(0)


    
    def intgrn(self):
        '''
		Creates a folder with all GRN interactions and how many times they appear
        in the expression files of the MultiDatasetNetworkAssembling object.
        
		Parameters:
		-----------
        None.
        
		Returns:
		--------
		None.
		
		Notes:
		------
		Creates a folder titled GRNX where X represents the number of files in the 
        MultiDatasetNetworkAssembling object. Inside the folder contains text files indicating
        how many times each interaction appeared and in what expression files it appeared in.
		
		Examples:
		---------
		>>> classInstance = MultiDatasetNetworkAssembling(numFiles = 4)
		>>> classInstance.intgrn()
		
        '''
        d = 'GRN_%iDatasets' % len(self.pls)
        try:
            os.mkdir(d)
        except:
            c = 1
            while ('%s_(%i)' % (d,c) in os.listdir('.')):
                c += 1
            d = '%s_(%i)' % (d,c)
            os.mkdir(d)
        data = {}
        edgepool = {}
        for i in range(len(self.pls)):
            pl = self.pls[i]
            print 'Reading %s...' % pl
            datatmp = {}
            text = open(self.grns[i]).read()
            for line in text.splitlines()[1:]:
                tmp = line.split('\t')
                edge = tmp[1] + ';;' + tmp[2]
                if edge not in edgepool.keys():
                    edgepool[edge] = 1
                else:
                    edgepool[edge] += 1
                datatmp[edge] = int(tmp[0])
            data[pl] = datatmp
        maxoverlap = 0
        overlaps = {}
        print 'Processing overlap...'

        for edge in edgepool:
            if edgepool[edge] > maxoverlap:
                maxoverlap = edgepool[edge]
            if edgepool[edge] not in overlaps.keys():
                overlaps[edgepool[edge]] = []
            overlaps[edgepool[edge]].append(edge)

        ranks = {}
        for o in overlaps.keys():
            print 'Processing overlap %i...' % o
            ranks = []
            for edge in overlaps[o]:
                r = []
                p = []
                for pl in self.pls:
                    if edge in data[pl]:
                        r.append(data[pl][edge])
                        p.append(pl)
                if len(r) != o:
                    raise KeyError
                rank = edge.split(';;')+[np.mean(r),str(p)]
                ranks.append(rank)
            ranks.sort(key=lambda ind:ind[2])
            text = ['Transcription Factor\tTarget\tAverage Rank\tPipelines']
            for rank in ranks:
                line = rank[0]
                for col in rank[1:]:
                    line += '\t%s' % str(col)
                text.append(line)
            open('%s/overlap_%i.txt' % (d,o),'w').write('\n'.join(text))
        self.igrn = d
    
    def intppi(self):
        '''
		Creates a folder with all PPI interactions and how many times and where they appear
        in the expression files of the MultiDatasetNetworkAssembling object.
        
		Parameters:
		-----------
        None.
        
		Returns:
		--------
		None.
		
		Notes:
		------
		Creates a folder titled PPIX where X represents the number of files in the 
        MultiDatasetNetworkAssembling object. Inside the folder contains text files indicating
        how many times each interaction appeared and in what expression files it appeared in.
		
		Examples:
		---------
		>>> classInstance = MultiDatasetNetworkAssembling(numFiles = 4)
		>>> classInstance.intppi()
		
        '''
        d = 'PPI_%iDatasets' % len(self.pls)
        try:
            os.mkdir(d)
        except:
            c = 1
            while ('%s_(%i)' % (d,c) in os.listdir('.')):
                c += 1
            d = '%s_(%i)' % (d,c)
            os.mkdir(d)
        data = {}
        edgepool = {}
        for i in range(len(self.pls)):
            pl = self.pls[i]
            print 'Reading %s...' % pl
            datatmp = {}
            text = open(self.ppis[i]).read()
            for line in text.splitlines()[1:]:
                tmp = line.split('\t')
                edge = tmp[1] + ';;' + tmp[2]
                if edge not in edgepool.keys():
                    edgepool[edge] = 1
                else:
                    edgepool[edge] += 1
                datatmp[edge] = int(tmp[0])
            data[pl] = datatmp
        maxoverlap = 0
        overlaps = {}
        print 'Processing overlap...'
        for edge in edgepool:
            if edgepool[edge] > maxoverlap:
                maxoverlap = edgepool[edge]
            if edgepool[edge] not in overlaps.keys():
                overlaps[edgepool[edge]] = []
            overlaps[edgepool[edge]].append(edge)

        ranks = {}
        for o in overlaps.keys():
            print 'Processing overlap %i...' % o
            ranks = []
            for edge in overlaps[o]:
                r = []
                p = []
                for pl in self.pls:
                    if edge in data[pl]:
                        r.append(data[pl][edge])
                        p.append(pl)
                if len(r) != o:
                    raise KeyError
                rank = edge.split(';;')+[np.mean(r),str(p)]
                ranks.append(rank)
            ranks.sort(key=lambda ind:ind[2])
            text = ['Transcription Factor\tTarget\tAverage Rank\tPipelines']
            for rank in ranks:
                line = rank[0]
                for col in rank[1:]:
                    line += '\t%s' % str(col)
                text.append(line)
            open('%s/overlap_%i.txt' % (d,o),'w').write('\n'.join(text))
        self.ippi = d
