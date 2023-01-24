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
cd ~

mkdir evan

cd evan

mkdir name_of_animal #here replace the name

cd name_of_animal
</font>

<font color='green'>jupyter notebook &</font>

In Jupyter Notebook open this ipynb program teplate and copy it with a different name.

## <font color='red'>2.Configuration and download section</font>

<font color='green'>species = 'cottoperca_gobio' # replace with your animal</font>

<font color='green'>!whereis pyensembl  # here you will see the path to species.py</font>

Edit **species.py** in directory /home/marta/anaconda3/envs/evan/lib/python3.11/site-packages/pyensembl/ #**change to your directory**

add in the end something like


<font color='green'>

anabas = Species.register(
    
    latin_name="anabas_testudineus",
    
    synonyms=["anabas"],
    
    reference_assemblies={
        
        "fAnaTes1.1": (94, MAX_ENSEMBL_RELEASE),
        
        "fAnaTes1.1": (108, MAX_ENSEMBL_RELEASE)
    })


cottoperca = Species.register(
    
    latin_name="cottoperca_gobio",
    
    synonyms=["cottoperca"],
    
    reference_assemblies={
        
        "fCotGob3.1": (108, MAX_ENSEMBL_RELEASE)
    })

</font>

### <font color='red'>Define the cache for downloaded files</font>

<font color='green'>!export PYENSEMBL_CACHE_DIR=~/.cache/pyensembl  # similar to this example </font>

**Define the Ensembl release you want to work with**

<font color='green'>!pyensembl install --release 108  --species cottoperca_gobio # replace the release and animal name </font>

Look at the .cache anf if not all required files are downloaded, you can do it manually with **wget**

<font color='green'>import wget</font>

<font color='green'>wget.download('https://ftp.ensembl.org/pub/release-108/fasta/cottoperca_gobio/dna/Cottoperca_gobio.fCotGob3.1.dna_sm.toplevel.fa.gz') # replace the path and animal name</font>

Verify if you are in the directory with the animal e.g. in this case "cottoperca_gobio"

<font color='green'>!pwd  # show the actual directory</font>

This file you must unzip it - see the example.

<font color='green'>!gzip -d /home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/Cottoperca_gobio/fCotGob3.1.dna_sm.toplevel.fa.gz  # replace path</font>

## <font color='red'>3.Import section</font>

<font color='green'>
import pyensembl

from Bio import SeqIO

from Bio.SeqUtils import gc_fraction

import matplotlib.pyplot as plt

import pandas as pd

import numpy as np
</font>

<font color='green'>files_108=!ls /home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/  # change the path of your own/font>

<font color='green'>files_108  # show files in your directory  </font>

Choose what file you will work with / e.g. <font color='green'>file_108[4]</font> if the file is 'Cottoperca_gobio.fCotGob3.1.dna_sm.toplevel.fa'

# DNA only relese ..... - pandas table

<font color='green'>window = 1000  # sliding windows - you can change the window size<font color='green'>

## GC DNA release .....

<font color='green'>

sec_values_108 = '' # replace with your release number

c=[]

for rec in SeqIO.parse(f"{where}"+files_108[4], "fasta"):
    
    c.append(rec.id)
    
    sec_values_108 += rec.seq
    
gc_values_108 = tuple(gc_fraction(sec_values_108[i:i+window]) for i in range(0,len(sec_values_108),window)) # replace with your 
release number

columns_dna_108 = ['GC DNA release 108'] # replace with your release number

gc_df_dna_108 = pd.DataFrame(gc_values_108, columns = columns_dna_108) # replace with your release number
</font>

## GC soft and (small) and unmasked (capital) repetitions

## GC DNA release .....

<font color='green'>

GC_scafolds_108={} # replace the release no.

for rec in SeqIO.parse("/home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/"+files_108[4], "fasta"):  # replace the path
    
    if not rec.id.startswith('CA'):  # here put the begin of the scafolds - to gain only chromosomes
        
        soft_mask=sum(1 for i in rec.seq if i in 'acgt')
        
        ALL_all=sum(1 for i in rec.seq if i in 'acgtACGT')
        
        GC_fraction=gc_fraction(rec.seq)
    
        GC_scafolds_108[rec.id]=(GC_fraction,soft_mask/ALL_all)  # we distract N and n   
</font>

<font color='green'>GC_scafolds_108 # replace with your release number</font>

You will see similar output

{'1': (0.40814745909188405, 0.17440378326633393),

 '2': (0.41938191102427813, 0.18024962740993744),
 
 '3': (0.4032972952361542, 0.1537007277720803),
 
 '4': (0.4064355195179001, 0.16729886597141377),
 
 '5': (0.40771028056836156, 0.18034140153997397),
 
 '6': (0.4092075763725821, 0.16868645101899657),
 
 '7': (0.40999156577005563, 0.15450076047047542),
 
 '8': (0.4115169054844562, 0.17200722014657274),
 
 '9': (0.40867709843459904, 0.1666861849883061),
 
 '10': (0.40778491354275276, 0.17446868825940207),
 
 '11': (0.4079062285611757, 0.16067698667209926),
 
 '12': (0.40598729214636164, 0.17566199710051028),
 
 '13': (0.40909998285626287, 0.1750210664566801),
 
 '14': (0.40489221966870464, 0.17078008966428682),
 
 '15': (0.409567413127844, 0.1616261429115718),
 
 '16': (0.40960814643715465, 0.1616937024780161),
 
 '17': (0.4080977423537412, 0.17158274561064024),
 
 '18': (0.41717795227880733, 0.16238750679625308),
 
 '19': (0.4116962929985457, 0.16647170628914018),
 
 '20': (0.41443687295086556, 0.1432771457382805),
 
 '21': (0.40558444142770966, 0.15142028127868343),
 
 '22': (0.41253309579159136, 0.16841521079242047),
 
 '23': (0.41754265112409306, 0.1597707466135346),
 
 '24': (0.41086758349515057, 0.18739326314781485)}

### Write the species data into csv file

<font color='green'>
import csv

with open('release_108_profile_soft_unmask_all_cottoperca.csv','w') as f:  # name csv file as you wish 
    
    w = csv.writer(f)
    
    w.writerows(GC_scafolds_108.items())
</font>

## Graph scatter - color by value

<font color='green'>

fig = plt.figure(figsize=(15,10))

ax1 = fig.add_subplot(111)

ax1.set_facecolor("lightgrey") 

#print(range(i,i+window),len(gc_values_94[i:i+window]))

        
names = list(GC_scafolds_108.keys())

v = list(GC_scafolds_108.values())

v1=[i[0] for i in v]

v2=[i[1] for i in v]

sc=ax1.scatter(names,v1, s=150, c=v2, cmap='RdYlGn', marker="o", ) # only GC fraction in y axe, only chromosomes

plt.grid(True)

plt.title('GC mean values of Cottoperca Gobio release 108')  # Replace the animal name

plt.ylabel('GC fraction')

plt.xlabel('Chromosomes')

plt.xticks(rotation = 45)

plt.colorbar(sc,label="Red 0% soft-masked to green 100% soft-masked", orientation="horizontal")

plt.savefig('release_108_profile_soft_unmask_per_chromosomes_cottoperca_gobio.png')  # Replace the animal name

plt.show()

</font>

## Relations between release ..... and .... - whole sequence

### GC DNA release 94

Replace the path with your

<font color='green'>

sec_values_94 = ''

c=[]

for rec in SeqIO.parse("/home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/"+files_108[4]"+files_94[4], "fasta"):  # replace the 
path

    sec_values_94 += rec.seq

gc_94_all = gc_fraction(sec_values_94)
</font>

### GC DNA release 108

<font color='green'>

sec_values_108 = ''

c=[]

for rec in SeqIO.parse("/home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/"+files_108[4], "fasta"): # replace the path

    sec_values_108 += rec.seq
    
gc_108_all = gc_fraction(sec_values_108)
</font>

## Results - %GC rel.94 / %GC rel.108 (not here)

## Results - here %GC rel.108 only

<font color='green'>
#print(f'%GC release 94 is {gc_94_all}, %GC release 108 is {gc_108_all}, %GC rel.94 / %GC rel 108 is {gc_94_all/gc_108_all}')

print(f'%GC release 108 is {gc_108_all}')
</font>

### Soft vs Unmasked

<font color='green'>

GC_all={}

'''

soft_mask=sum(1 for i in sec_values_94 if i in 'acgt')

ALL_all=sum(1 for i in sec_values_94 if i in 'acgtACGT')

GC_fraction=gc_fraction(sec_values_94)

GC_all['94']=(GC_fraction,soft_mask/ALL_all)  # distract N and n   

'''

soft_mask=sum(1 for i in sec_values_108 if i in 'acgt')

ALL_all=sum(1 for i in sec_values_108 if i in 'acgtACGT')

GC_fraction=gc_fraction(sec_values_108)

GC_all['108']=(GC_fraction,soft_mask/ALL_all)  # distract N and n , replace release no.
</font>

<font color='green'>
GC_all.items() # show result
</font>

### GC to ATGC, gc to atgc

comment what is not appropriate / here rel.94

<font color='green'>

GC_to_ATGC={}

'''

atgc_lower=sum(1 for i in sec_values_94 if i in 'atgc')

gc_lower=sum(1 for i in sec_values_94 if i in 'gc')

atgc_upper=sum(1 for i in sec_values_94 if i in 'ATGC')

gc_upper=sum(1 for i in sec_values_94 if i in 'GC')

GC_to_ATGC['94_lower']=gc_lower/atgc_lower  # distract N and n   

GC_to_ATGC['94_upper']=gc_upper/atgc_upper  # distract N and n 

'''

atgc_lower=sum(1 for i in sec_values_108 if i in 'atgc')

gc_lower=sum(1 for i in sec_values_108 if i in 'gc')

atgc_upper=sum(1 for i in sec_values_108 if i in 'ATGC')

gc_upper=sum(1 for i in sec_values_108 if i in 'GC')

GC_to_ATGC['108_lower']=gc_lower/atgc_lower  # distract N and n, replace release no.

GC_to_ATGC['108_upper']=gc_upper/atgc_upper  # distract N and n, replace release no.
</font>

<font color='green'>GC_to_ATGC.items()</font>

<font color='green'>
You will see 

dict_items([('108_lower', 0.4053924169972052), ('108_upper', 0.4105225443714316)])
</font>

## GC soft and (small) and unmasked (capital) repetitions - by chromosomes

## GC DNA release ..... - by chromosomes

<font color='green'>import numpy as np</font>

<font color='green'>window = 1000 # choose the window size</font>

<font color='green'>files_108[4]  # Check if you choose right file</font>

<font color='green'>
for rec in SeqIO.parse("/home/marta/.cache/pyensembl/fCotGob3.1/ensembl108/"+files_108[4], "fasta"): # replace the path

    GC_scafolds_108_windowed={}  # per chromosome
    
    if not rec.id.startswith('CA'):
    
        #print(rec.id)  
        
        for w in range(0,len(rec.seq),window):
        
              
            #print(range(i,i+window),len(gc_values_94[i:i+window]))
            
            soft_mask=sum(1 for i in rec.seq[w:w+window] if i in 'acgt')
            
            ALL_all=sum(1 for i in rec.seq[w:w+window] if i in 'acgtACGT')
            
            GC_fraction=gc_fraction(rec.seq[w:w+window])
            
            try:
            
                GC_scafolds_108_windowed[rec.id+'|'+str(w)]=(GC_fraction,soft_mask/ALL_all)  # distract N and n 
                
                #print(range(i,i+window),len(gc_values_94[i:i+window])) 
                
            except ZeroDivisionError:
            
                GC_scafolds_108_windowed[rec.id+'|'+str(w)]=(GC_fraction,0.0)  # distract N and n
            
        fig = plt.figure(figsize=(50,20))
        
        ax1 = fig.add_subplot(111)
        
        ax1.set_facecolor("lightgrey")   
        
        names = list(GC_scafolds_108_windowed.keys())
        
        v = list(GC_scafolds_108_windowed.values())
        
        v1=[i[0] for i in v]
        
        v2=[i[1] for i in v]
            
        sc=ax1.scatter(names,v1, s=5, c=v2, cmap='RdYlGn', marker="o", ) # only GC fraction in y axe, only chromosomes
        
        ax1.xaxis.set_ticks(np.arange(0, len(rec.seq)//1000, len(rec.seq)//12000),rotation = 45)
        
        ax1.grid(True)        

        plt.title(f'GC mean values of Cottoperca Gobio release 108 - chromosome {rec.id}') # replace the animal name and release no.
        plt.ylabel('GC fraction')
        
        plt.xlabel(f'Chromosome {rec.id} windows')
        
        #plt.xticks(rotation = 45)
        
        plt.colorbar(sc,label="Red 0% soft-masked to green 100% soft-masked", orientation="horizontal")
        
        plt.savefig(f'release_108_profile_soft_unmask_{rec.id}.png')  # replace the animal name and release no.
        
        plt.show()
        
        
        with open(f'release_108_profile_soft_unmask_{rec.id}.csv','w') as f: # replace the animal name and release no.
        
            csv_file = csv.writer(f)
            
            csv_file.writerows(GC_scafolds_108_windowed.items())
                 
</font>        

# Put all png images to pdf file

In the end all the files are put together into the pdf file.

<font color='green'>
import os

from fpdf import FPDF
</font>

<font color='green'>imagelist=!ls | grep .png</font>
<font color='green'>
del imagelist[-1] # if needy, delete what does not delong to chromosomes, here the last item from the list
</font><font color='green'>
imagelist=sorted(imagelist,key= lambda x:int(x.split('_')[-1].split('.')[0])) # sort the chromosomes acc. to their numbers
</font>
<font color='green'>
pdf = FPDF(orientation = 'L')

# imagelist is the list with all image filenames

for image in imagelist:

    pdf.add_page()
    pdf.image(image,x=1,y=10,w=300,h=135)
    
pdf.output("Release108_Cottoperca_gobio_graphs.pdf", "F")  # replace the name and relese no.
</font>


```python

```
