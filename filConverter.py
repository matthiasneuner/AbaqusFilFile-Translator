import sys
import os
import numpy as np
import math

import utils.exportEngine as eE
from utils.inputfileparser import parseInputFile, printKeywords

if __name__ == "__main__":
    
    if len(sys.argv) < 3 or sys.argv[1] == '--help':
        print('Usage: filConverter.py  FILFILE.fil  EXPORTDEFINITION.inp')
        print('')
        print('Available Keywords:')
        print('')
        printKeywords()
        
    fn = sys.argv[1]
    jobFile = sys.argv[2]
    
    exportJobs = parseInputFile(jobFile)
    
    exportName  = ''.join(fn.split('/')[-1].split('.')[-2])
    print(fn)
    print("{:<20}{:>20}".format('opening file',fn))
    print('*'*40)
    wordsize = 8  
    
    exportEngine = eE.ExportEngine(exportJobs, exportName )
    
    chunkSize = 513* wordsize
    batchSize = chunkSize * 4096 * 32 # = ~ 538 MByte
    fileStat = os.stat(fn)
    fileSize = fileStat.st_size
    
    numberOfBatchSteps = math.ceil(fileSize / batchSize)
    
    print("file has a size of {:} bytes".format(fileSize))
    print("file will be processed in {:} steps".format(numberOfBatchSteps))
    
    batchIdx = 0
    currentIndex = 0
    fn = np.memmap(fn, dtype='b', mode='r+', )
    while batchIdx < fileSize:
        fileRemainder = fileSize - batchIdx
        idxEnd = batchIdx + (batchSize if fileRemainder >= batchSize else fileRemainder)
        fil = np.copy(fn[batchIdx :   idxEnd])
        words = fil.reshape( -1 , chunkSize )
        words = words[:, 4:-4]
        words = words.reshape(-1, 8)
        
        while currentIndex < len(words):
            recordLength = eE.filInt(words[currentIndex])[0]  
            if recordLength<=2:
                print('found a record with 0 length content, possibly aborted Abaqus analysis')
                break
            if currentIndex + recordLength > len(words):
                batchIdx += int(math.floor(currentIndex/512))* 513 * 8  # move to beginning of the current 512 word block in the batchChunk and restart with a new bathChunk
                currentIndex =  ( (currentIndex%512) )                  # of course, restart at the present index
                break
            recordType = eE.filInt(words[currentIndex+1])[0]
            recordContent = words[currentIndex+2 : currentIndex+recordLength]
            success = exportEngine.computeRecord(recordLength, recordType, recordContent)
            currentIndex += recordLength
        if currentIndex == len(words):
            currentIndex = 0
            batchIdx += batchSize
        del fil, words
            
    exportEngine.finalize()
        
    print('*'*40)
    print('Summary of .fil file:')
    print('{:<20}{:>20}'.format('nodes:',len(exportEngine.abqNodes)))
    print('{:<20}{:>20}'.format('elements:',len(exportEngine.abqElements)))
    print('{:<20}{:>20}'.format('element sets:',len(exportEngine.abqElSets)))
    for setName, elList in exportEngine.abqElSets.items():
        print('{:.<4}{:<16}{:>11} elements'.format('.',setName, len(elList)))
    print('{:<20}{:>20}'.format('increments:',exportEngine.nIncrements))

    

