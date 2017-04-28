import sys
import os
import numpy as np
import math

import utils.exportEngine as eE
from utils.inputfileparser import parseInputFile, printKeywords
import time


if __name__ == "__main__":
    
    if len(sys.argv) < 3 or sys.argv[1] == '--help':
        print('Usage: filConverter.py  FILFILE.fil  EXPORTDEFINITION.inp')
        print('')
        print('Available Keywords:')
        print('')
        printKeywords()
        exit()
        
    fn = sys.argv[1]
    jobFile = sys.argv[2]
    
    exportJobs = parseInputFile(jobFile)
    
    exportName  = ''.join(fn.split('/')[-1].split('.')[-2])
    lockFile =  fn.split('.')[0]+'.lck' 
    print("{:<20}{:>20}".format('opening file',fn))
    print('*'*48)
    wordsize = 8  
    
    exportEngine = eE.ExportEngine(exportJobs, exportName )
    
    chunkSize = 513* wordsize
    batchSize = chunkSize * 4096 * 32  # = ~ 538 MByte  ... size in BYTES
    fileStat = os.stat(fn)
    fileSize = fileStat.st_size
    
    numberOfBatchSteps = math.ceil(fileSize / batchSize)
    
    print("file has a size of {:} bytes".format(fileSize))
    print("file will be processed in {:} steps".format(numberOfBatchSteps))
    
    fileIdx = 0
    wordIdx = 0

    parseFile = True
    while parseFile:
        try:
            fileStat = os.stat(fn)
            fileSize = fileStat.st_size
            
            if fileIdx < fileSize:
                fnMap = np.memmap(fn, dtype='b', mode='r', )
                fileRemainder = fileSize - fileIdx # remaining file size in BYTES
                idxEnd = fileIdx + (batchSize if fileRemainder >= batchSize else fileRemainder) #get end index 
                
                # in case we are operating on an unfinished file and 'catch' an unfinished chunk
                idxEnd -= idxEnd % chunkSize
                                   
                batchChunk = np.copy(fnMap[fileIdx :   idxEnd])  # get chunk of file
                words = batchChunk.reshape( -1 , chunkSize )  # get words 
                words = words[:, 4:-4]
                words = words.reshape(-1, 8)
                while wordIdx < len(words):
                    recordLength = eE.filInt(words[wordIdx])[0]  
                    if recordLength<=2:
                        print('found a record with 0 length content, possible an aborted Abaqus analysis')
                        if os.path.exists(lockFile):
                            print("found .lck file, waiting for new result .fil data")
                            time.sleep(5)
                            fnMap = np.memmap(fn, dtype='b', mode='r', )
                            batchChunk = np.copy(fnMap[fileIdx :   idxEnd])  # get chunk of file
                            words = batchChunk.reshape( -1 , chunkSize )  # get words 
                            words = words[:, 4:-4]
                            words = words.reshape(-1, 8)
                            continue
                        else:
                            parseFile = False
                            break
                    
                    # the next record exceeds our batchChunk, so we do some trick:
                        # - set the wordIdx to the end of the so far progressed frame
                        # - move the frame to the wordIDx
                    if wordIdx + recordLength > len(words):
                        bytesProgressedInCurrentBatch = int(math.floor(wordIdx/512))* 513 * 8
                                                           
                                                           
                        if bytesProgressedInCurrentBatch == 0: # indicator for an aborted analysis
                            print('terminated file, possible an aborted Abaqus analysis')
                            if os.path.exists( lockFile):
                                print("found .lck file, waiting for new result .fil data")
                                time.sleep(5)
                                fnMap = np.memmap(fn, dtype='b', mode='r', )
                                batchChunk = np.copy(fnMap[fileIdx :   idxEnd])  # get chunk of file
                                words = batchChunk.reshape( -1 , chunkSize )  # get words 
                                words = words[:, 4:-4]
                                words = words.reshape(-1, 8)
                                continue
                            else:
                                parseFile = False
                                break
                        
                        fileIdx +=   bytesProgressedInCurrentBatch# move to beginning of the current 512 word block in the batchChunk and restart with a new bathChunk
                        wordIdx =  ( (wordIdx%512) )                  # of course, restart at the present index
                        break
                    
                    recordType = eE.filInt(words[wordIdx+1])[0]
                    recordContent = words[wordIdx+2 : wordIdx+recordLength]
                    success = exportEngine.computeRecord(recordLength, recordType, recordContent)
                    wordIdx += recordLength
                    
                # clean finish of a batchChunk
                if wordIdx == len(words):
                    wordIdx = 0
                    fileIdx = idxEnd
#                    fileIdx += batchSize
                    
                del words
            else:
                if os.path.exists( lockFile ):
                    print("found .lck file, waiting for new result .fil data or CTRL-C to finish...")
                    time.sleep(10)
                else:
                    break
            
        except KeyboardInterrupt:
            print("Interrupted by user")
            break
            
    exportEngine.finalize()
        
    print('*'*48)
    print('Summary of .fil file:')
    print('{:<20}{:>28}'.format('nodes:',len(exportEngine.allNodes)))
    print('{:<20}{:>28}'.format('elements:',len(exportEngine.allElements)))
    print('{:<20}{:>28}'.format('element sets:',len(exportEngine.elSets)))
    for setName, elSet in exportEngine.elSets.items():
        for elType, elements in elSet.elements.items():
            print('{:.<4}{:<16}{:8}{:>11} elements'.format('.',setName,elType, len(elements)))
    print('{:<20}{:>28}'.format('increments:',exportEngine.nIncrements))

    

