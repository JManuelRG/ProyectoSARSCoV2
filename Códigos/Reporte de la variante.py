# -*- coding: utf-8 -*-
"""
@author: José Manuel
"""

# importar todos lo métodos y librerías requeridos
from dna import dna 
import numpy as np 
import matplotlib.pyplot as plt 
from scov import numpy_image_dict
from helper import *

# Leer secuencia de la variante
# Cambiar nombre de la variante
nombre='B.1.1.159_Mx.txt'
dict_seq_1 = read_dna_seq(nombre)

genes=[]
longitudes=[]
# Generar reporte
print(f"Longitud de la secuencia por proteína del genoma {nombre}")
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
           
print(f"Longitud total = {sum(ordenlong)}")
# Gráfico del reporte
#https://datatofish.com/horizontal-bar-chart-matplotlib/
# Barras horizontales:
fig, ax = plt.subplots()
bars = ax.barh(ordengenes, ordenlong)
ax.bar_label(bars)    
ax.barh(ordengenes,ordenlong, color='green')
ax.set_title(f'Longitud por Gen de la secuencia {nombre}')
ax.set_ylabel('Genes')
ax.set_xlabel('Longitud')
ax.text(17500,'ORF10',s=f"Total={sum(ordenlong)}")
plt.show()
fig.savefig(fname=f"Reporte de secuencia {nombre[:-5]}.png", 
            dpi=1024)

