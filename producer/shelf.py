import shelve

YukawaDict = {}

for line in open('data/mhmodp_13000_higgs_A_highmA.txt', 'r'):
    elements = line.rsplit()

    tanb = float(elements[0])
    mass = float(elements[1])
    Yt = float(elements[5])
    Yb = float(elements[6])

    if mass == 1050:        
        savestr = 'tanb_' + str(tanb).replace('.0', '')
        YukawaDict[savestr] = {'Yt':Yt, 'Yb':Yb}


db = shelve.open('Yukawa_A.db')
db['Yukawa'] = YukawaDict
db.close()
