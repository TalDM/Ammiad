import pandas
import glob, os
import re
from itertools import izip

path = ''
path2 = ""
outputdir='divided/'

myFastqs = pandas.read_excel(path+"1992_plate.xlsx")
fnames=list();
R2_fnames=list();
trimGC = 4

for columns,j in myFastqs.iterrows():
    filename=list()
    currentSerial=j['Serial']
    currentlib=j['library']
    #currentplant=j['plant']
    finalplant=j['finalplant']
    currentsublib=j['sub_library']
    currentBarcode=j['Barcode']
    gene=j['gene']
    currenttrim=j['trim']


    treatID=str(finalplant)
    print (treatID)

    for filename in glob.glob(os.path.join(path,'*R1*.fastq')):
        R2_filename=re.sub("R1","R2",filename)

        if filename.startswith(path2+'Undetermined'):
            print("nice!")
            print(filename)
            print(R2_filename)
            fnames.append(filename)
            R2_fnames.append(R2_filename)
    print ("working on these files: "+str(fnames))
    print ("looking for this barcode: "+currentBarcode)
    newname=outputdir+treatID+"_R1.fastq"
    newnameR2=outputdir+treatID+"_R2.fastq"
    f= open(newname, 'w')
    f_R2= open(newnameR2, 'w')

    count = 0

    if gene == 'GBS':
        print (gene)
        print (currentlib)
        print (currenttrim)

    for i in range(len(fnames)):
        name=fnames[i]
        nameR2=R2_fnames[i]
        fhand = open(name)
        R2_fhand = open(nameR2)

    with open(name) as fhand:

        for line1, line1_R2 in izip(fhand, R2_fhand):

            read_identifiers=line1
            read_seq=fhand.next()
            plusLine=fhand.next()
            read_qual=fhand.next()

            read_identifiers_R2=line1_R2
            read_seq_R2=R2_fhand.next()
            plusLine_R2=R2_fhand.next()
            read_qual_R2=R2_fhand.next()

            temp=read_seq.startswith(currentBarcode)
            if read_seq.startswith(currentBarcode):

                count+=1
                f.write(read_identifiers+read_seq[currenttrim:]+plusLine+read_qual[currenttrim:])
                f_R2.write(read_identifiers_R2+read_seq_R2[trimGC:]+plusLine_R2+read_qual_R2[trimGC:])

    f.close()
    f_R2.close()
    fnames=list();
    R2_fnames=list()
