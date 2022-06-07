# importar los modulos requeridos
from dna import dna
import numpy as np
import matplotlib.pyplot as plt
from scov import numpy_image_dict
from helper import *

# Virus Variant
variant='B.1.1.159_Mx.txt'

# Read the dna sequence file-1 previously downloaded from NCBI.
dict_seq_1 = read_dna_seq(variant)
# Modify the sequence with dummy 'N' nucleotide.
dict_seq_1 = gene_mod(dict_seq_1)

# Read the dna sequence file-2 previously downloaded from NCBI.
dict_seq_2 = read_dna_seq('Wuhan_2019Dec.txt')
# Modify the sequence with dummy 'N' nucleotide.
dict_seq_2 = gene_mod(dict_seq_2)

# Variables para determinar longitudes de los genes en orden
#Longitudes en orden    
genes=[]
longitudes=[]

for key in dict_seq_1:
    genes.append(dict_seq_1[key][0])
    longitudes.append(len(dict_seq_1[key][1]))
    print(dict_seq_1[key][0],":",len(dict_seq_1[key][1]))    

ordengenes=['ORF1ab',
            'S',
            'ORF3a',
            'E',
            'M',
            'ORF6',
            'ORF7a',
            'ORF7b',
            'ORF8',
            'N',
            'ORF10']

#Longitudes en orden    
ordenlong=[len(dict_seq_1['gene=ORF1ab'][1]),
           len(dict_seq_1['gene=S'][1]),
           len(dict_seq_1['gene=ORF3a'][1]),
           len(dict_seq_1['gene=E'][1]),
           len(dict_seq_1['gene=M'][1]),
           len(dict_seq_1['gene=ORF6'][1]),
           len(dict_seq_1['gene=ORF7a'][1]),
           len(dict_seq_1['gene=ORF7b'][1]),
           len(dict_seq_1['gene=ORF8'][1]),
           len(dict_seq_1['gene=N'][1]),
           len(dict_seq_1['gene=ORF10'][1])]


# Definir función para listar el aminoácido:
def num2nucleot(i):
    var=""
    if i == 0:
        var='A'
    if i == 255:
        var='T'
    if i == 100:
        var='C'
    if i == 200:
        var='G'
    if i == 75:
        var='N'
    return var 

# Diccionario
mutaciones={}

    
# Create matplotlib subplots for each gene. 
f,ax = plt.subplots(nrows=11,ncols=3,figsize=(25,30))
gene_name = list(numpy_image_dict.keys())
row = 0
col = 0
mut_dict={}
for i in gene_name:
    G = i[5:]
    # Loop thru each gene in the Cornona Virus nucleotide sequence.
    gene_us = dna(dict_seq_1['gene='+G][1])
    # Invoke the transcription method of the class dna 
    gene_us.transcription()
    # Invoke the mothod that converts the gene sequence into a numpy array.
    numpfy_usa = gene_us.numpfy()
    # Reshape the numpy array with a predeifned shape from the numpy_image_dict dictionary.
    numpfy_usa = numpfy_usa.reshape(numpy_image_dict['gene='+G][0])
    # sub-plot the numpy array with matplotlib pcolor method.
    ax[row][col].pcolor(numpfy_usa)
    ax[row][col].set_title(G+' '+variant[:-4])
    col+=1
    gene_china = dna(dict_seq_2['gene='+G][1])
    # Invoke the transcription method of the class dna 
    gene_china.transcription()
    # Invoke the mothod that converts the gene sequence into a numpy array.
    numpfy_china = gene_china.numpfy()
    # Reshape the numpy array with a predeifned shape from the numpy_image_dict dictionary.
    numpfy_china = numpfy_china.reshape(numpy_image_dict['gene='+G][0])
    # sub-plot the numpy array with matplotlib pcolor method.
    ax[row][col].pcolor(numpfy_china)
    ax[row][col].set_title(G+' Wuhan')
    col+=1

    # To find the gene mutation subtract the numpy array from base sequence with the newer sequence. Here the 
    # the Chinese sequence is the base sequence and the USA sequence is a newer sequence.
    mut = numpfy_china - numpfy_usa
    if mut.any():
        # Here we are looking for a non zero value in the mutated numpy array (result of the subtracting the 2 numpy arrays).
        # Presence of non-zero value means that there is difference between the 2 numpy arrays and the gene has 
        # mutataions. If there are mutations in the gene create a python dictionary "mut_dict" with details as below.
        # {'<Gene_Name-1>': [[<value_of_base_seq>, <value_of_newer_seq>, <value_in_mutated_numpy>, (x_value,y_value)]], '<Gene_Name-2>': [[<value_of_base_seq>, <value_of_newer_seq>, <value_in_mutated_numpy>, (x_value,y_value)]]}
        mut_nec = np.nonzero(mut)
        x=mut_nec[0]
        y=mut_nec[1]
        l=0
        mut_dict[G]=[]
        for i in x:
            us_base = numpfy_usa[i][y[l]]
            ch_base = numpfy_china[i][y[l]]
            mut_base = mut[i][y[l]]
 
            #Agregado
            nucl_us=num2nucleot(us_base) #Transforma el número al nucleótido
            nucl_ch=num2nucleot(ch_base)
            info_list = [us_base,ch_base,mut_base,(i,y[l])]
            mut_dict[G].append(info_list)
            #Sustituido por el reporte final de abajo
            #print("Nucleotido {} en Wuhan mutado a {} en Mx_B.1.1 en la posicion {} en el Gen {}".format(nucl_us,nucl_ch,(i,y[l]),G))
            #print("Nucleotido {} en Wuhan mutado a {} en {} en la posicion {} en el Gen {}".format(nucl_us,nucl_ch,variant[:-4],((i+1)*(y[l]+1)),G))
            nucl_us,nucl_ch,((i+1)*(y[l]+1)),G
            
            l+= 1
    # Giving a title to the matplotlib subplot
    ax[row][col].pcolor(mut)
    ax[row][col].set_title(G+' Mutaciones')
    row+= 1
    col=0

f.tight_layout()
# Saving the matplotlib subplot as a png.
f.savefig(f'Sars_Cov-2_Mutaciones en {variant[:-4]}.png')


### Generar Reporte Final ordenado de mutaciones
#info_list= [us_base,ch_base,mut_base,(i,y[l])]
mut_report=dict()
mutorden=dict() # Almacena las mutaciones a reportar
for i in mut_dict.keys(): # Iterar en cada gen
    mut_report[i]=[] # Crear las claves para cada gen vacías
    mutorden[i]=[]
    for fila in range(len(mut_dict[i])): #Iterar en cada lista
        mut_report[i].append([num2nucleot(mut_dict[i][fila][0]),
                           num2nucleot(mut_dict[i][fila][1]),
                           (mut_dict[i][fila][3][0]+1)*(mut_dict[i][fila][3][1]+1)])
        mutorden[i].append(mut_report[i][fila][2])

# Ordenar reporte de mutaciones
for i in mutorden.keys():
    mutorden[i].sort()


print("Reporte de mutaciones:")
# Imprimir reporte final:        
for i in mut_report.keys(): # Iterar en cada gen
    print(f"{len(mut_report[i])} Mutacion(es) detectada(s) en el gen {i}")
    print(f"Longitud total del gen {i}={ordenlong[ordengenes.index(i)]}")
    for fila in range(len(mut_report[i])): #Iterar en cada lista
        print(f"     La base {mut_report[i][fila][1]} (Wuhan) mutó a {mut_report[i][fila][0]} en {variant} en la posición {mut_report[fila][2]}")
        
