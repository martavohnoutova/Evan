# Profiling of GC of species in Ensembl releases

<font color='red'>**Guide of usage**</font>

date: 31.12.2022

Marta Vohnoutová

version 0.2

## <font color='red'>1. Introduction</font>

The goal of this tool is to profile GC% and repeats proportions along entire chromosome sequences.

The concrete configuration steps follow.

### <font color='red'>1.1 Setup of environment</font>

This is the template to run your own setup of GC profiling of your species.

a. **At first, you need to create your environment.**

- **install Anaconda (description for linux Ubuntu v.22)**

<font color= 'green'>

$ sudo apt update

$ sudo apt install curl

$ curl --output anaconda.sh https://repo.anaconda.com/archive/Anaconda3-5.3.1-Linux-x86_64.sh

$ sha256sum Anaconda3-5.3.1-Linux-x86_64.sh

$ bash Anaconda3-5.3.1-Linux-x86_64.sh


</font>

b. **Activate your environment**

<font color= 'green'>conda create -n evan biopython</font>

<font color='green'>source activate evan</font>

c. **Install necessary packages**

<font color='green'>conda install scipy matplotlib jupyter-notebook pip pandas numpy wget pyensembl gzip fpdf</font>

d. **Start Jupyter Notebook**

Go to the home directory, and start Juputer Notebook.

<font color='green'>


## <font color='red'>2.Import section</font>

<font color='green'>

#import pyensembl

from Bio import SeqIO

from Bio.SeqUtils import gc_fraction

import matplotlib.pyplot as plt

import pandas as pd

import numpy as np

import os

import csv

from collections import OrderedDict

from PIL import Image

from glob import glob

</font>

## <font color='red'>3.Configuration and download section </font>

### <font color='red'> 3.1 Variable part - replace with right values </font>

<i>
species = 'cottoperca_OIST' # replace with your animal

where_i_am = '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/' # replace with your path

where_is_fasta = '/home/rsymonova/Evan/Evan/isochores/' #replace with the path to fasta files

window = 1000  # sliding windows - you can change the window size

column_name = 'GC DNA' # replace with your column name

what_we_filter = 'not rec.id.isdigit()'

</i>


<i>
filesRS=!ls {where_is_fasta} | grep .fa  #replace grep with right filter

filesRS

['Anabas_OIST.fa',

 'Cottoperca_OIST.fa',
 
 'Danio.fa',
 
 'Esox.fa',
 
 'Lepisosteus.fa',
 
 'Petromyzon.fa',
 
 'Salmo.fa',
 
 'Xiphophorus.fa']
 
my_fasta_file=filesRS[1]  # replace with your choice
</i>

Choose what file you will work with / e.g. file_108[4] if the file is 'Cottoperca_gobio.fCotGob3.1.dna_sm.toplevel.fa'

### <font color='red'>Create subdirectory for particular species</font>
<i>
try:

    os.mkdir(f'{where_i_am}{species}')
    
except OSError as error:

    pass
</i>

#where are you

<i>
!pwd

/mnt/data/radka/Evan/Vohnoutova_et_al2023
</i>

If the file you must unzip it - see the example.

!gzip -d /home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/Cottoperca_gobio/fCotGob3.1.dna_sm.toplevel.fa.gz  # replace the path

where = '/home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/'

</i>

## <font color='red'>3.Parsing section
​
### <font color='red'>DNA only relese ..... - pandas table</i>

### <font color='red'>GC DNA release .....

sec_values = '' # replace with your release number

<i>
c=[]

for rec in SeqIO.parse(f"{where_is_fasta}{my_fasta_file}", "fasta"):
    
    c.append(rec.id)
    
    sec_values += rec.seq
    
gc_values = tuple(gc_fraction(sec_values[i:i+window]) for i in range(0,len(sec_values),window)) # replace with your release number
​
columns_dna = [column_name] # replace with your name e.g. release number
​
gc_df_dna = pd.DataFrame(gc_values, columns = columns_dna) # replace with your release number

gc_values

</i>

### <font color='red'>GC soft and (small) and unmasked (capital) repetitions</font>

### <font color='red'>GC DNA release .....</font>

<i>

GC_chromosoms=OrderedDict()

chomosome_no=0

for rec in SeqIO.parse(f"{where_is_fasta}{my_fasta_file}", "fasta"):  # replace the path

    if eval(what_we_filter):  # chromosomes only
    
        soft_mask=sum(1 for i in rec.seq if i in 'acgt')
        
        ALL_all=sum(1 for i in rec.seq if i in 'acgtACGT')
        
        GC_fraction=gc_fraction(rec.seq)
        
        chomosome_no += 1
        
        GC_chromosoms[str(chomosome_no) + '+' + rec.id]=((GC_fraction,soft_mask/ALL_all),len(rec.seq))  # we distract N and n   
​
#all_images = sorted([os.path.basename(f) for f in glob(f"{where_i_am}{species}/*profile*.png")],key=lambda x: int(x.split('_')[-1].split('.')[0]))
​
GC_chromosoms = OrderedDict(sorted(GC_chromosoms.items(), key=lambda x: x[1][1],reverse=True))

GC_chromosoms # replace with your release number

OrderedDict([('5+ENA|LR131935|LR131935.1',

              ((0.40771028056836156, 0.26172963543221006), 30479438)),
             ('9+ENA|LR131939|LR131939.1',
             
              ((0.40867709843459904, 0.23503259514623812), 30072547)),
             ('3+ENA|LR131933|LR131933.1',
             
              ((0.4032972952361542, 0.24405482959290087), 30028995)),
             ('4+ENA|LR131934|LR131934.1',
             
              ((0.4064355195179001, 0.26340883704278506), 28949045)),
             ('13+ENA|LR131920|LR131920.1',
             
              ((0.40909998285626287, 0.2553741620635248), 27742916)),
             ('6+ENA|LR131936|LR131936.1',
             
              ((0.4092075763725821, 0.2419895700504481), 27680345)),
             ('10+ENA|LR131917|LR131917.1',
             
              ((0.40778491354275276, 0.2757633966163528), 27438269)),
             ('1+ENA|LR131916|LR131916.1',
             
              ((0.40814745909188405, 0.27957934285544367), 27055867)),
             ('16+ENA|LR131923|LR131923.1',
             
              ((0.40960814643715465, 0.23370477984345692), 26581398)),
             ('14+ENA|LR131921|LR131921.1',
             
              ((0.40489221966870464, 0.250162566794175), 25704503)),
             ('17+ENA|LR131924|LR131924.1',
             
              ((0.4080977423537412, 0.2534343737321996), 25156145)),
             ('15+ENA|LR131922|LR131922.1',
             
              ((0.409567413127844, 0.23324127625295674), 24963210)),
             ('21+ENA|LR131929|LR131929.1',
             
              ((0.40558444142770966, 0.24334879786129016), 24104793)),
             ('8+ENA|LR131938|LR131938.1',
              ((0.4115169054844562, 0.26902241688608086), 23432536)),
             ('7+ENA|LR131937|LR131937.1',
             
              ((0.40999156577005563, 0.2396669825521754), 23069680)),
             ('12+ENA|LR131919|LR131919.1',
             
              ((0.40598729214636164, 0.2718863070962306), 22902521)),
             ('22+ENA|LR131930|LR131930.1',
             
              ((0.41253309579159136, 0.252774641384023), 22609112)),
             ('24+ENA|LR131932|LR131932.1',
             
              ((0.41086758349515057, 0.2701858233286377), 22441820)),
             ('11+ENA|LR131918|LR131918.1',
             
              ((0.4079062285611757, 0.2485296108056439), 22193419)),
             ('19+ENA|LR131926|LR131926.1',
             
              ((0.4116962929985457, 0.25594310984477486), 21060126)),
             ('20+ENA|LR131928|LR131928.1',
             
              ((0.41443687295086556, 0.2043045254502776), 17599799)),
             ('23+ENA|LR131931|LR131931.1',
             
              ((0.41754265112409306, 0.22188361473468693), 15932363)),
             ('18+ENA|LR131925|LR131925.1',
             
              ((0.41717795227880733, 0.22128985444522287), 14930383)),
             ('2+ENA|LR131927|LR131927.1',
             
              ((0.41938191102427813, 0.23513998467310082), 12923129))])
              
with open(f'{where_i_am}{species}/{species}.csv','w') as f:  # name csv file as you wish 

    w = csv.writer(f)
    
    w.writerows(GC_chromosoms.items())
</i>

### <font color='red'>Graph scatter - color by value</i>
<i>
fig = plt.figure(figsize=(15,10))

ax1 = fig.add_subplot(111)

ax1.set_facecolor("lightgrey")     

#print(range(i,i+window),len(gc_values_94[i:i+window]))
        
names = list(i.split('+')[0] for i in GC_chromosoms.keys())

v = list(GC_chromosoms.values())

v1=[i[0][0] for i in v]

v2=[i[0][1] for i in v]
​
sc=ax1.scatter(names,v1, s=150, c=v2, cmap='RdYlGn', marker="o", ) # only GC fraction in y axe, only chromosomes

plt.grid(True)

plt.title(f'GC mean values of {species}\nsorted by size')  # Replace the animal name

plt.ylabel('GC fraction')

plt.xlabel('Chromosomes')

plt.xticks(rotation = 45)

plt.colorbar(sc,label="Red 0% soft-masked to green 100% soft-masked", orientation="horizontal")

plt.savefig(f'{where_i_am}{species}/soft_unmask_per_chromosomes_{species}.png')  # Replace the animal name

plt.show()
</I>
​
Relations between release ..... and .... - whole sequence
GC DNA release 94
Replace the path with your

sec_values_94 = ''
c=[]
for rec in SeqIO.parse("/home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/"+files_108[4]"+files_94[4], "fasta"):  # replace the path
    sec_values_94 += rec.seq
​
gc_94_all = gc_fraction(sec_values_94)
​
GC DNA release 108
sec_values = ''
c=[]
for rec in SeqIO.parse(f"{where_is_fasta}{my_fasta_file}","fasta"): # replace the path
    sec_values += rec.seq
    
gc_all = gc_fraction(sec_values)
Results - %GC rel.94 / %GC rel.108 (not here)
Results - here %GC rel.108 only
#print(f'%GC release 94 is {gc_94_all}, %GC release 108 is {gc_108_all}, %GC rel.94 / %GC rel 108 is {gc_94_all/gc_108_all}')
print(f'%GC of {species} is {gc_all}')
%GC of cottoperca_OIST is 0.4092285307602729
Soft vs Unmasked
gc_all.items()
dict_items([('cottoperca_OIST', (0.4092285307602729, 0.25015370660549935))])
GC to ATGC, gc to atgc
comment what is not appropriate / here rel.94

GC_to_ATGC={}
​
atgc_lower=sum(1 for i in sec_values if i in 'atgc')
gc_lower=sum(1 for i in sec_values if i in 'gc')
atgc_upper=sum(1 for i in sec_values if i in 'ATGC')
gc_upper=sum(1 for i in sec_values if i in 'GC')
​
GC_to_ATGC['lower']=gc_lower/atgc_lower  # distract N and n, replace release no. 
GC_to_ATGC['upper']=gc_upper/atgc_upper  # distract N and n, replace release no.
​
GC_to_ATGC.items()
dict_items([('lower', 0.4023479510007686), ('upper', 0.4115239378931707)])
GC soft and (small) and unmasked (capital) repetitions - by chromosomes
GC DNA release ..... - by chromosomes
my_fasta_file  # Check if you choose right file
'Cottoperca_OIST.fa'
lengths=[]
for rec in SeqIO.parse(f"{where_is_fasta}{my_fasta_file}", "fasta"): 
        lengths.append(len(rec.seq))
max_len=max(lengths)
print(max_len)
30479438
max_len//12000
2539
# adjusted chromosome density
for rec in SeqIO.parse(f"{where_is_fasta}{my_fasta_file}", "fasta"): # replace the path
    GC_chromosoms_windowed={}  # per chromosome
    if eval(what_we_filter):
        #print(rec.id)  
        for w in range(0,len(rec.seq),window):
              
            #print(range(i,i+window),len(gc_values_94[i:i+window]))
            soft_mask=sum(1 for i in rec.seq[w:w+window] if i in 'acgt')
            ALL_all=sum(1 for i in rec.seq[w:w+window] if i in 'acgtACGT')
            GC_fraction=gc_fraction(rec.seq[w:w+window])
            
            try:
                GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,soft_mask/ALL_all)  # distract N and n    
                #print(range(i,i+window),len(gc_values_94[i:i+window])) 
            except ZeroDivisionError:
                GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,0.0)  # distract N and n
            
        #fig = plt.figure(figsize=(int(50* (len(rec.seq)/max_len)),int(20* (len(rec.seq)/max_len))))
        fig = plt.figure(figsize=(50,20))
        ax1 = fig.add_subplot(111)
        ax1.set_facecolor("lightgrey")   
        
        names = list(GC_chromosoms_windowed.keys())
        
        v = list(GC_chromosoms_windowed.values())
        v1=[i[0] for i in v]
        v2=[i[1] for i in v]
        names_part=[n.split('|')[0] for n in names]
       
        sc=ax1.scatter(names_part,v1, s=5, c=v2, cmap='RdYlGn', marker="o", ) # only GC fraction in y axe, only chromosomes
        ax1.xaxis.set_ticks(np.arange(0, max_len//1000, max_len//12000))
        ax1.grid(True)        
​
        plt.title(f'GC% values of {species} - chromosome {rec.id}', fontsize = 35) # replace the animal name and release no.
        plt.ylabel('GC fraction', fontsize = 30)
        plt.xlabel(f'Chromosome {rec.id} windows', fontsize = 30)
        plt.xticks(fontsize = 30, rotation = 45)
        plt.yticks(fontsize = 30)
        plt.colorbar(sc,label="Red 0% soft-masked to green 100% soft-masked", orientation="horizontal")
        plt.savefig(f'{where_i_am}{species}/{species}_profile_soft_unmask_{rec.id}.png')  # replace the animal name and release no.
        plt.show()
        
        
        with open(f'{where_i_am}{species}/{species}_profile_soft_unmask_{rec.id}.csv','w') as f: # replace the animal name and release no.
            csv_file = csv.writer(f)
            csv_file.writerows(GC_chromosoms_windowed.items())
                 
        
























Merge all graph images
Prepare small graphs without colorbars
# adjusted chromosome density
for rec in SeqIO.parse(f"{where_is_fasta}{my_fasta_file}", "fasta"): # replace the path
    GC_chromosoms_windowed={}  # per chromosome
    if eval(what_we_filter):
        #print(rec.id)  
        for w in range(0,len(rec.seq),window):
              
            #print(range(i,i+window),len(gc_values_94[i:i+window]))
            soft_mask=sum(1 for i in rec.seq[w:w+window] if i in 'acgt')
            ALL_all=sum(1 for i in rec.seq[w:w+window] if i in 'acgtACGT')
            GC_fraction=gc_fraction(rec.seq[w:w+window])
            
            try:
                GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,soft_mask/ALL_all)  # distract N and n    
                #print(range(i,i+window),len(gc_values_94[i:i+window])) 
            except ZeroDivisionError:
                GC_chromosoms_windowed[str(w)+'|'+ rec.id]=(GC_fraction,0.0)  # distract N and n
            
        #fig = plt.figure(figsize=(int(50* (len(rec.seq)/max_len)),int(20* (len(rec.seq)/max_len))))
        fig = plt.figure(figsize=(25,10))
        ax1 = fig.add_subplot(111)
        
        names = list(GC_chromosoms_windowed.keys())
        v = list(GC_chromosoms_windowed.values())
        v1=[i[0] for i in v]
        v2=[i[1] for i in v]
        names_part=[n.split('|')[0] for n in names]
       
        sc=ax1.scatter(names_part,v1, s=5, c=v2, cmap='RdYlGn', marker="o", ) # only GC fraction in y axe, only chromosomes
        ax1.xaxis.set_ticks(np.arange(0, max_len//1000, max_len//12000))    
        plt.title(f'GC% values of {species} - chromosome {rec.id}', fontsize = 35) # replace the animal name and release no.        
        plt.savefig(f'{where_i_am}{species}/{species}_nobar_soft_unmask_{rec.id}.{len(rec.seq)}.png')  # replace the animal name and release no.
        plt.show()
























Read graphs one by one and merge them all
# list of all graph images sorted to chromosomes no.
all_images = sorted([f for f in glob(f"{where_i_am}{species}/*nobar*.png")],key=lambda x: x.split('.')[-2], reverse = True)
all_images
['/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131933|LR131933.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131916|LR131916.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131920|LR131920.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131924|LR131924.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131935|LR131935.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131921|LR131921.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131936|LR131936.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131923|LR131923.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131930|LR131930.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131934|LR131934.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131939|LR131939.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131938|LR131938.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131917|LR131917.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131932|LR131932.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131918|LR131918.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131929|LR131929.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131937|LR131937.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131926|LR131926.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131919|LR131919.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131922|LR131922.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131928|LR131928.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131931|LR131931.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131925|LR131925.1.png',
 '/home/rsymonova/Evan/Evan/Vohnoutova_et_al2023/cottoperca_OIST/cottoperca_OIST_nobar_soft_unmask_ENA|LR131927|LR131927.1.png']
images = [Image.open(f"{x}") for x in all_images]
widths, heights = zip(*(i.size for i in images))
max_width = max(widths)
total_height = sum(heights)
new_im = Image.new('RGB', (max_width, total_height))
​
y_offset = 0
for im in images:
  new_im.paste(im, (0,y_offset,))
  y_offset += im.size[1]
​
new_im.save(f'{where_i_am}{species}/{species}_all_graphs.png')
Read small graphs one by one and merge them all to two columns
def get_concat(all_images,im_width,im_height):
    
    images = [Image.open(all_images[i]) for i in range(len(all_images))]
    if len(all_images)%2 == 0:
        how_high=len(all_images)
    else:
        how_high=len(all_images)+1
        
    dst = Image.new('RGB', (2*im_width, int((how_high * im_height)/2)))
    for i in range(0,len(images)+1,2):
        try:
            dst.paste(images[i], (0, int((i * im_height)/2)))
        except IndexError:
            pass
        try:            
            dst.paste(images[i+1], (im_width,  int((i * im_height)/2)))
        except IndexError:
            pass
    return dst
im = Image.open(all_images[0])
im_width=im.size[0]
im_height=im.size[1]
im.close()
get_concat(all_images,im_width,im_height).save(f'{where_i_am}{species}/{species}_all_small_graphs.png')
Put all png images to pdf file
In the end all the files are put together into the pdf file.

import os
from fpdf import FPDF
imagelist=!ls | grep .png
del imagelist[-1] # delete what does not delong to chromosomes
imagelist=sorted(imagelist,key= lambda x:int(x.split('_')[-1].split('.')[0]))
pdf = FPDF(orientation = 'L')
# imagelist is the list with all image filenames
for image in imagelist:
    pdf.add_page()
    pdf.image(image,x=1,y=10,w=300,h=135)
pdf.output("Profiles_gastropod.pdf", "F")  # replace the name and relese no.
​

