# Author: Alejandro V. Seba
# Creation Date:
## Description:

'''
Checks the risk haplotype in the haplotype freq files
coming from plink output in script haplotypes_pleiotropy.sh
'''
import sys

def readFreqHapFile(path):
    '''
    '''
    hapFreqs={}
    with open(path, 'r') as fhtp:
        for line in fhtp:
            (key, val) = line.split(' ')
            hapFreqs[(key)] = float(val.strip())
    return hapFreqs

def checkRiskHaploFreqs(hapFreqs, riskHap, pleio, agon_early_early, agon_late_late, agon_early_late,
antagon_early_early, antagon_late_late, antagon_early_late, OnsetA, OnsetB, threeshold):
    '''
    '''
    if len(hapFreqs) == 2:
#        print '2 Haplotypes'
        if riskHap in hapFreqs:
            classify_onset_disease(OnsetA, OnsetB, threeshold, agon_early_early,
             agon_late_late, agon_early_late, pleio)
        else:
            classify_onset_disease(OnsetA, OnsetB, threeshold, antagon_early_early,
             antagon_late_late, antagon_early_late, pleio)
    else:
#        print len(hapFreqs), 'Haplotypes'
        justFreqs = sorted(hapFreqs.values())[:-2]
        minim = sum(justFreqs)
        if minim > 0.1:
            #print(hapFreqs)
            sys.exit(['This is exceeding the frequency threshold.'])
        elif riskHap in hapFreqs and hapFreqs[riskHap] > 0.1:
#            print 'AGON'
            classify_onset_disease(OnsetA, OnsetB, threeshold, agon_early_early,
             agon_late_late, agon_early_late, pleio)
        else:
#            print 'ANT'
            classify_onset_disease(OnsetA, OnsetB, threeshold, antagon_early_early,
             antagon_late_late, antagon_early_late, pleio)

#He eliminado el filtro de edad debido a que las infecciosas
#tienen una edad inferior a la esperada (1 año por defecto)
def writeFiles(agon_early_early, agon_late_late, agon_early_late,
antagon_early_early, antagon_late_late, antagon_early_late, path_out):
    '''
    '''
    print(path_out)
    if agon_early_early:
        with open(path_out + "Agon_early_early", 'a') as agF:
            for i in agon_early_early:
                agF.write(str(i))
                agF.write('\n')
    if agon_late_late:
        with open(path_out + "Agon_late_late", 'a') as agF:
            for i in agon_late_late:
                agF.write(str(i))
                agF.write('\n')
    if agon_early_late:
        with open(path_out + "Agon_early_late", 'a') as agF:
            for i in agon_early_late:
                agF.write(str(i))
                agF.write('\n')
    if antagon_early_early:
        with open(path_out + "Antagon_early_early", 'a') as agF:
            for i in antagon_early_early:
                agF.write(str(i))
                agF.write('\n')
    if antagon_late_late:
        with open(path_out + "Antagon_late_late", 'a') as agF:
            for i in antagon_late_late:
                agF.write(str(i))
                agF.write('\n')
    if antagon_early_late:
        with open(path_out + "Antagon_early_late", 'a') as agF:
            for i in antagon_early_late:
                agF.write(str(i))
                agF.write('\n')


def classify_onset_disease(OnsetA, OnsetB, threeshold, list1, list2, list3, pleio):
    if float(OnsetA) <= threeshold and float(OnsetB) <= threeshold:
        list1.append(pleio)
    elif float(OnsetA) > threeshold and float(OnsetB) > threeshold:
        list2.append(pleio)
    elif (float(OnsetA) <= threeshold and float(OnsetB) > threeshold) or \
    (float(OnsetA) > threeshold and float(OnsetB) <= threeshold):
        list3.append(pleio)


def main():
    """
    main function
    """
    agon_early_early = []
    agon_late_late = []
    agon_early_late = []
    antagon_early_early = []
    antagon_late_late = []
    antagon_early_late = []
#     path = '/home/jrodriguez/Projects/PLEIOTROPY/Hapls_Output/rs4766578_rs3184504.fhtp'
    path = sys.argv[1]
    age_threeshold = int(sys.argv[2])
#    riskHap='AC' ## ${riskHap}
    riskHap = sys.argv[3]
    print(riskHap)
#     car='rs1333042	Coronary heart disease	A	rs6475606	Intracranial aneurysm	T	LATE	EARLY'
    OnsetA = sys.argv[8]
    OnsetB = sys.argv[13]
    pleio_rel = '\t'.join(sys.argv[4:14]) #pleiotropic relation from haplotypes_pleiotropies.sh ${CAR}

    out_path = sys.argv[-1]
## Input of pleiotropic relation from haplotypes_pleiotropies.sh ${CAR}
    #age = sys.argv[4] #age threshold
# LD Threshold/file analyzed name
    #ld = sys.argv[5]
    hapFreqs = readFreqHapFile(path)

    checkRiskHaploFreqs(hapFreqs, riskHap, pleio_rel, agon_early_early, agon_late_late, agon_early_late,
    antagon_early_early, antagon_late_late, antagon_early_late, OnsetA, OnsetB, age_threeshold)

    writeFiles(agon_early_early, agon_late_late, agon_early_late,
    antagon_early_early, antagon_late_late, antagon_early_late,
    out_path)

if __name__ == "__main__":
    exit(main())
