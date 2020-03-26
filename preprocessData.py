import numpy as np
import subprocess
import json
import os
from os.path import isfile, join
from collections import defaultdict
import matplotlib.pyplot as plt
import scipy 
from scipy.spatial.distance import euclidean 

def getSNPs(settings, path):
    # Get the SNPS for each case/control file as specified
    # by "path" and then finds the common SNPs (intersection).
    # Output: list of common SNPs to cases or controls files.

    # Get .bim files of controls
    files = [file for file in os.listdir(path) if isfile(join(path, file)) and '.bim' in file]

    # Initialize
    totalSNPs = list()
    SNPs = list()

    # Geet SNPs
    for file in files:
        print("Parsing SNPs from " + file)
        with open(join(path, file), 'r') as inFile:
            for row in inFile:
                row = row.split('\t')
                SNPs.append(row[1])
            totalSNPs.append(SNPs)
    print('Finished parsing')

    # Get intersection
    commonSNPs = set(totalSNPs[0])
    for i in range(1, len(totalSNPs)):
        commonSNPs = commonSNPs.intersection(set(totalSNPs[i]))
    commonSNPs = list(commonSNPs)

    return commonSNPs

def filterSNPs(settings):
    # Filters the SNPs for each case/control file and 
    # only considers the SNPs common to all files (cases and controls)
    # Output: binary files copied to FiltBinaries

    print('Filtering common SNPs')
    subprocess.call(['bash', 'filterSNPs.sh'])
    print('Common SNPs filtered')

def getMergelist(settings):
    # Creates a merge file to be used as input for plink 
    # to merge files (cases and controls)

    # Get control files
    filePath = settings['directory']['GWAS_controlsFILTBIN']
    files = [f for f in os.listdir(filePath) if isfile(join(filePath,f)) and '.bim' in f]
    files = [f.split('.')[0] for f in files]

    # Write mergelist
    with open(settings['file']['mergeList.txt'],'w') as outFile:
        for f in files:
            outFile.write(f + '.bed ' + f + '.bim ' + f + '.fam' + '\n')


def mergeFiles(settings):
    # Uses the mergeList to merge files and creates 
    # a file called LGI1 in Data/GWAS

    # Merge files
    print('Merging files')
    subprocess.call(['bash', 'mergeFiles.sh'])
    print('Files successfully merged')

def IBDfilt(settings):
    # This function computes and IBD computation and creates 
    # A list of subjects to remove removeIDs.txt in Resources 

    # Run IBD computation
    if bool(settings['IBD_perform']):
        print('Compute IBD')
        subprocess.call(['bash','IBD.sh'])
        print("IBD successfully computed")
    else:
        print('IDB already computed')

    # Get list of patients to be removed 
    IBD_IDs = list()
    with open(settings['file']['IBD'], 'r') as inFile:
        next(inFile)
        for row in inFile:
            row = row.split('\t')
            row = [r for r in row if r is not '']
            PI_hat = row[8]
            if float(PI_hat) > settings['IBD_threshold']:
                IBD_IDs.append(row[0])

    # Input genotypes
    f = settings['file']['LGI1']
    sep = '\t'
    addPhenotype(settings, f, sep)           

    # Get list of patients and cases 
    dataset = list()
    with open('./Data/modLGI1.fam','r') as inFile:
        for row in inFile:
            row = row.split('\t')
            if int(row[5].split('\n')[0]) != -9:
                dataset.append(row[0])
    
    # Get list of IBDS patient/cases
    high_IBD = [ID for ID in np.unique(IBD_IDs) if ID in dataset]

    # Load controls and cases
    controls = list() ; cases = list()
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile: controls.append(row.split('\"')[1])
    with open(settings['file']['GWASIDsCases'], 'r') as inFile: 
        for row in inFile: cases.append(row.split('\"')[1])
    
    # Check number of cases and controls
    case_count = 0 ; control_count = 0 
    for ID in high_IBD:
        if ID.split('_')[-1].split('.')[0] in cases:
            case_count +=1
        elif ID.split('_')[-1].split('.')[0] in controls:
            control_count += 1
    print('Number of cases to be excluded: ' + str(case_count))
    print('Number of controls to be excluded: ' + str(control_count))

    # Create exclusion file 
    with open(settings['file']['excludeID.txt'], 'w') as outFile:
        for ID in high_IBD: 
            outFile.write(ID + ' ' + ID + ' ' + '\n')


def LGI1QC(settings):
    # Computes quality control in the merged files 
    # 1. Creates temporary files with the rsIDs correpsonding to MHC CHR 6 
    # 2. Filters out the MHC CHR 5 and maf at 0.05 and geno and missingness at 0.1
    # Output: filtLGI1 in ./Data

    # Perform quality control
    print("Performing quality control")
    subprocess.call(['bash', 'qualityControl.sh'])
    print("Quality control successfully performed")

def LGI1PCA(settings):
    # Computes PCA on LGI1 master files with common SNPs betwen cases and controsl 
    # no MHC and maf at 0.05 and geno and missingness at 0.1

    if bool(settings['performPCA']) is True:
        # Perform PCA
        print("Performing PCA")
        subprocess.call(['bash', 'pca.sh'])
        print("PCA successfully performed")
    else:
        print("PCA deactivated")

def plotPCA(settings, f):

    # Initialize
    pc1 = list()
    pc2 = list()
    badge = list()

    # Read file
    with open(f, 'r') as inFile:
        for row in inFile:
            row = row.split(' ')
            pc1.append(float(row[2]))
            pc2.append(float(row[3]))
            badge.append(int(row[1].split('_')[1]))

    # Plot
    pc1 = np.array(pc1) ; pc2 = np.array(pc2) ; badge =  np.array(badge)
    fig, ax = plt.subplots()
    col = list()
    count = 0 
    for b in np.unique(badge):
        idx = np.where(badge == b)
        col.append(np.random.rand(3,))
        ax.scatter(x = pc1[idx], y = pc2[idx], c = col[count], \
            label = b, alpha = 0.3, linewidths=1)   
        count += 1
    plt.title("Principal components")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend()
    plt.show(block=False)

    # Plot PCs individually 

    fig1, ax1 = plt.subplots(nrows=3, ncols=7)
    fig1.subplots_adjust(hspace=0.5)
    count = 0
    for ax, b in zip(ax1.flatten(), np.unique(badge)):
        idx = np.where(badge == b)
        ax.scatter(x = pc1[idx], y = pc2[idx], c = col[count], \
            label = b, alpha = 1, linewidths = 1)
        ax.set(title = "Plate " + str(b), xlabel = 'PC1', ylabel = 'PC2')
        count += 1
    plt.show()

def addPhenotype(settings, f, sep):
    # Opens filtLGI1 and adds the phenotype based on the list of 
    # controls and cases. If not in cases or controls, ID is added 
    # to list and written in exclude.txt
    # Output modLGI1.fam 

    # Read GWAS IDS
    controlsID = list() 
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile:
            controlsID.append(row.split('"')[1])
    casesID = list()
    with open(settings['file']['GWASIDsCases'], 'r') as inFile:
        for row in inFile:
            casesID.append(row.split('"')[1])

    # Feed phenotype to fam file and create exclude list
    exclude = list()
    with open(settings['file']['modLGI1.fam'], 'w') as outFile:
        with open(f + '.fam', 'r') as inFile:
            for row in inFile:
                row = row.split(sep)
                subjectID = row[0].split('_')[3].split('.')[0]
                if subjectID in controlsID:
                    row[-1] = '1'
                    outFile.write('\t'.join(row) + '\n')
                elif subjectID in casesID:
                    row[-1] = '2'
                    outFile.write('\t'.join(row) + '\n')
                else:
                    exclude.append(row[0])
                    outFile.write('\t'.join(row))

def refiltLGI1(settings):
    # Takes original data and filters as filtLGI1 but including MHC
    # Moves files to Data

    print("Refiltering and preparing data")
    subprocess.run(['bash','refilter.sh'])
    print("Refiltering finished")

    # Add phenotype again 
    f = settings['file']['filtLGI1']
    sep = ' '
    addPhenotype(settings, f, sep)


def patientMatching(settings):

    # Perform patient matching
    print("Computing patient matching as the Euclidean distance between each case and control")

    # GET IDS and PCs
    patIDs = np.array([])
    patIDsGWAS = np.array([])
    PCs = list()
    IDs = list()
    with open(settings['file']['PCAeigenvectors'], 'r') as inFile:
        for row in inFile:
            patIDsGWAS = np.append(patIDsGWAS, row.split(' ')[0])
            ID = np.asanyarray(row.split(' ')[0].split('_')[-1].split('.')[0])
            IDs.append(row.split(' ')[0].split('_')[-1].split('.')[0])
            patIDs = np.append(patIDs, ID)
            PCs.append(np.asanyarray([float(pc) for pc in row.split(' ')[2:]]))

    # Open cases ID
    casesID = np.array([])
    with open(settings['file']['GWASIDsCases'], 'r') as inFile:
        for row in inFile:
            casesID = np.append(casesID, row.split("\"")[1])

    # Open controls ID and get the ones in PatIDs (some removed in filtering)
    controlsID = np.array([])
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile:
            controlsID = np.append(controlsID, row.split("\"")[1])
    controlsID = np.asanyarray([id for id in controlsID if id in patIDs])

    # Get matching patients
    ratio = settings['ControlCaseRatio']

    # For each case, compute euclidean distance
    keepIndxs = list()
    for case in casesID:

        # If case filtered or in exclusion list, continue 
        if case not in patIDs:
            continue
        
        # Get the index and PCs of the case. Store.
        idx = int(np.where(patIDs == case)[0][0])
        keepIndxs.append(idx)
        patPCs = PCs[idx]

        # Compute euclidean distance against all controls 
        patEuDist = list()
        for control in controlsID: 

            # If control filtered or in exclusion list, continue
            if control not in patIDs:
                continue

            # Get controls PCs and compute euclidean distance
            ctrlidx = int(np.where(patIDs == control)[0][0])
            ctrlPCs = PCs[ctrlidx]
            patEuDist.append(euclidean(patPCs, ctrlPCs))

        # Get lower distances between case and controls
        tempidx = np.asanyarray(patEuDist).argsort()[:ratio]
        ctrlsIDtemp = controlsID[tempidx]

        # Keep indxs
        for ctrl in ctrlsIDtemp:
            keepIndxs.append(np.where(patIDs == ctrl)[0][0])
        
    keepIndxs = np.unique(keepIndxs)
    # Write list of patients 
    patIDsGWAS = patIDsGWAS[keepIndxs]
    with open(settings['file']['patList'], 'w') as outFile:
        for pat in patIDsGWAS:
            outFile.write(pat + ' ' + pat + '\n')

def postMatchPCA(settings):

    print("Performing PCA after matching")
    subprocess.call(['bash','postPCA.sh'])
    print("PCA finished")
    f = settings['file']['postPCAeigenvectors']
    plotpostPCA(settings, f)

def plotpostPCA(settings, f):

     # Initialize
    pc1 = list()
    pc2 = list()
    IDs = list()

    # Read file
    with open(f, 'r') as inFile:
        for row in inFile:
            row = row.split(' ')
            pc1.append(float(row[2]))
            pc2.append(float(row[3]))
            IDs.append(str(row[1].split('_')[-1].split('.')[0]))
    
    # Load cases and controls ID
    controlsID = list() 
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile:
            controlsID.append(row.split('"')[1])
    casesID = list()
    with open(settings['file']['GWASIDsCases'], 'r') as inFile:
        for row in inFile:
            casesID.append(row.split('"')[1])

    # Get label (case vs control)
    label = list()
    for ID in IDs:
        if ID in casesID:
            label.append('cases')
        elif ID in controlsID:
            label.append('controls')

    # Plot
    pc1 = np.array(pc1) ; pc2 = np.array(pc2) ; label =  np.array(label)
    fig, ax = plt.subplots()
    col = list()
    count = 0 
    for b in np.unique(label):
        idx = np.where(label == b)
        col.append(np.random.rand(3,))
        ax.scatter(x = pc1[idx], y = pc2[idx], c = col[count], \
            label = b, alpha = 0.7, linewidths=1)   
        count += 1
    plt.title("Principal components")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend()
    plt.show(block=False)

    # Plot PCs individually 

    fig1, ax1 = plt.subplots(nrows=1, ncols=2)
    fig1.subplots_adjust(hspace=0.5)
    count = 0
    for ax, b in zip(ax1.flatten(), np.unique(label)):
        idx = np.where(label == b)
        ax.scatter(x = pc1[idx], y = pc2[idx], c = col[count], \
            label = b, alpha = 1, linewidths = 1)
        ax.set(title = "Plate " + str(b), xlabel = 'PC1', ylabel = 'PC2')
        count += 1
    plt.show()

def filtPostMatchData(settings):

    # Initialize
    PCs = list()
    IDs = list()
    f = settings['file']['postPCAeigenvectors']

    # Import PCs
    with open(f, 'r') as inFile:
        for row in inFile:
            row = row.split(' ')
            PCs.append([float(pc) for pc in row[2:]])
            IDs.append(row[1:2])

    # Filter PCs using thresholds
    idxs = list()
    count = 0
    for PC in PCs:
        if PC[0] > -0.00349 and PC[1] > 0.0554:
            count += 1
            continue
        elif PC[0] > 0.02377 and PC[0] < -0.075:
            count += 1
            continue
        elif PC[0] > 0.125:
            count += 1
            continue
        else:
            idxs.append(count)
            count += 1

    # Rewrite PCA 
    PCs = list(list(PC) for PC in np.asanyarray(PCs)[idxs]); 
    IDs = list(ID[0] for ID in np.asanyarray(IDs)[idxs])
    with open(settings['file']['filtposPCAeigenvectors'], 'w') as outFile:
        for PC, ID in zip(PCs, IDs):
            PC = [str(x) for x in PC]
            outFile.write(ID + ' ' + ID + ' ' + ' '.join(PC) + '\n')

    # Get list of patients and write 
    with open(settings['file']['filtPostMatchPatList'], 'w') as outFile:
        for pat in IDs:
            outFile.write(pat + ' ' + pat + '\n')

    # Plot PCA 
    f = settings['file']['filtposPCAeigenvectors']
    plotpostPCA(settings, f)


def associatioAnalysis(settings):

    # Perform association analysis
    print("Performing Asssociation analysis")
    subprocess.call(['bash','associationAnalysis.sh'])
    print ('\n')
    print("Association analysis performed successfully")


def main():

    # Open settings
    with open('settings.json', 'r') as inFile:
        settings = json.load(inFile)

    # # Convert to binary if not converted yet
    # if bool(settings['convert2Binary']) is True:
    #     print('Converting files to binary')
    #     subprocess.run(['bash', 'convert2Binary.sh'], stdout=subprocess.PIPE())
    # else:
    #     print("Files already converted to binary")

    # # Get common SNPs within controls and cases
    # controls_SNPs = getSNPs(settings, settings['directory']['GWAS_controlsBIN'])
    # cases_SNPs = getSNPs(settings, settings['directory']['GWAS_casesBIN'])

    # # Get common SNPs between cases and controls
    # commonSNPs = list(set(controls_SNPs).intersection(set(cases_SNPs)))

    # # Write common SNPS
    # with open(settings['file']['commonSNPs.txt'], 'w') as outFile:
    #     for snp in commonSNPs:
    #         outFile.write(snp + '\n')

    # # # Filter snps
    # filterSNPs(settings)
   
    # # Create merge list
    # getMergelist(settings)

    # # Merge files
    # mergeFiles(settings)

    # # Compute IBD and create list of patients to exclude 
    # IBDfilt(settings)

    # # Compute quality control
    # LGI1QC(settings)

    # # Add phenotype
    # f = settings['file']['filtLGI1']
    # sep = ' '
    # addPhenotype(settings, f, sep)

    # # Compute PCA to study badge biases
    # LGI1PCA(settings)

    # # Plot PCA
    # if bool(settings['plotPCA']):
    #     f = settings['file']['PCAeigenvectors']
    #     plotPCA(settings, f)

    # # Refilter to create master file 
    # refiltLGI1(settings)

    # # Patient matching
    # patientMatching(settings)

    # # Plot PCA
    # postMatchPCA(settings)

    # # Manually remove postMatch data
    # filtPostMatchData(settings)

    # # Add field separator in PCA
    # with open(settings['file']['PCAeigenvectors'], 'r') as inFile: 
    #     with open('./Resources/PCAtst.eigenvec', 'w') as outFile:
    #         header = ['PC' + str(x) for x in range(1,21)]
    #         header = ['FID','IID'] + header
    #         outFile.write('\t'.join(header) + '\n') 
    #         for row in inFile:
    #             row = row.split(' ')
    #             row = '\t'.join(row)
    #             outFile.write(row)

    # # Perform association analysis
    # associatioAnalysis(settings)


if __name__ == "__main__":
    main()    # Plot PCA
    
