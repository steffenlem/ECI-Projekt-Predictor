#contains different endcodings for peptide sequences

allAminoAcids = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']

blomap = {
    'A': [-0.57, 0.39, -0.96, -0.61, -0.69],
    'R': [-0.40, -0.83, -0.61, 1.26, -0.28],
    'N': [-0.70, -0.63, -1.47, 1.02, 1.06],
    'D': [-1.62, -0.52, -0.67, 1.02, 1.47],
    'C': [0.07, 2.04, 0.65, -1.13, -0.39],
    'Q': [-0.05, -1.50, -0.67, 0.49, 0.21],
    'E': [-0.64, -1.59, -0.39, 0.69, 1.04],
    'G': [-0.90, 0.87, -0.36, 1.08, 1.95],
    'H': [0.73, -0.67, -0.42, 1.13, 0.99 ],
    'I': [0.59, 0.79, 1.44, -1.90, -0.93 ],
    'L': [0.65, 0.84, 1.25, -0.99, -1.90],
    'K': [-0.64, -1.19, -0.65, 0.68, -0.13],
    'M': [0.76, 0.05, 0.06, -0.62, -1.59],
    'F': [1.87, 1.04, 1.28, -0.61, -0.16],
    'P': [-1.82, -0.63, 0.32, 0.03, 0.68],
    'S': [-0.39, -0.27, -1.51, -0.25, 0.31],
    'T': [-0.04, -0.30, -0.82, -1.02, -0.04],
    'W': [1.38, 1.69, 1.91, 1.07, -0.05],
    'Y': [1.75, 0.11, 0.65, 0.21, -0.41],
    'V': [-0.02, 0.30, 0.97, -1.55, -1.16]
}
test = [['A','A'],['V','V']]

def twentybit(thearray):
    translatedArray = []
    for x in thearray:
        innerArray = []
        for y in x:
            singleAA = [0]*20
            singleAA[allAminoAcids.index(y.upper())] = 1
            innerArray.extend(singleAA)
        translatedArray.append(innerArray)
    return translatedArray

def blomapCode(thearray):
    translatedArray = []
    for x in thearray:
        innerArray = []
        for y in x:
            innerArray.extend(blomap[y.upper()])
        translatedArray.append(innerArray)
    return translatedArray



print(twentybit(test))
print(blomapCode(test))
