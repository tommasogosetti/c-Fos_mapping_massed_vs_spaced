# -*- coding: utf-8 -*-
import os

def get_list_of_file_in_directory(dir): #ottiene lista file in una directory escludendo le altre directory ivi contenute (e in generale tutto ciò che non è file)
    if(not os.path.isdir(dir)):
        raise Exception(f"{dir} is not a directory")
    files = []
    for f in os.scandir(dir):
        if(f.is_file()):
            files.append(f.path)
    return files

file_list=get_list_of_file_in_directory('274_strutture')

with open('abbreviazioni.txt', 'r') as f:
    abbreviazioni=f.readlines()

for i in file_list:
    lista_righe_nuove=[]
    with open(i, 'r') as f:
        righe=f.readlines()
    with open(os.getcwd()+'\con_sigle\SIGLE_'+i.split('\\')[1], 'a') as f:
        for l in range(0,274):
            riga_splittata=righe[l].split(';')
            riga_splittata[0]=abbreviazioni[l].rstrip('\n')
            nuova_riga=';'.join(riga_splittata)
            f.write(nuova_riga)
        
print(os.getcwd())
