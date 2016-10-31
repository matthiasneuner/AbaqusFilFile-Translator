import sys
import numpy as np

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
    
    exportName  = ''.join(fn.split('.')[:-1])
    print("{:<20}{:>20}".format('opening file',fn))
    print('*'*40)
    wordsize = 8  
    
    exportEngine = eE.ExportEngine(exportJobs, exportName )
    
    with open(fn, 'rb') as fh:
        fil = np.fromfile(fn, dtype='b')
        words = fil.reshape( -1 , 513*8 )
        words = words[:, 4:-4]
        words = words.reshape(-1, 8)
        
        currentIndex = 0
        
        while currentIndex < len(words):
            recordLength = eE.filInt(words[currentIndex])[0]                
            recordType = eE.filInt(words[currentIndex+1])[0]
            recordContent = words[currentIndex+2 : currentIndex+recordLength]
            success = exportEngine.computeRecord(recordLength, recordType, recordContent)
            currentIndex += recordLength
        exportEngine.finalize()
        
    print('*'*40)
    print('Summary of .fil file:')
    print('{:<20}{:>20}'.format('nodes:',len(exportEngine.abqNodes)))
    print('{:<20}{:>20}'.format('elements:',len(exportEngine.abqElements)))
    print('{:<20}{:>20}'.format('element sets:',len(exportEngine.abqElSets)))
    for setName, elList in exportEngine.abqElSets.items():
        print('{:.<4}{:<16}{:>11} elements'.format('.',setName, len(elList)))
    print('{:<20}{:>20}'.format('increments:',exportEngine.nIncrements))

    

