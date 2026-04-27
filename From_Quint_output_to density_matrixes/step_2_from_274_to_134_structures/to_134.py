# -*- coding: utf-8 -*-
import statistics as stat
import os


def get_list_of_file_in_directory(dir): #ottiene lista file in una directory escludendo le altre directory ivi contenute (e in generale tutto ciò che non è file)
    if(not os.path.isdir(dir)):
        raise Exception(f"{dir} is not a directory")
    files = []
    for f in os.scandir(dir):
        if(f.is_file()):
            files.append(f.path)
    return files

def rename(lista_righe, first_line, last_line, sigla, nome_completo):
    riga=lista_righe[first_line].split(';')
    riga[0]=sigla
    riga[1]=nome_completo
    riga=';'.join(riga)
    print(riga)
    return riga
    
def unite(lista_righe, first_line, last_line, sigla, nome_completo):
    tot_objects=0
    tot_area=0
    avrg_density=0
    stdevs_list=[]
    avrg_stdev=0
    for i in range(first_line, last_line+1):
        sottoarea_splittata=lista_righe[i].split(';')
        tot_objects+=float(sottoarea_splittata[2])
        tot_area+=float(sottoarea_splittata[3])
        stdevs_list.append(float(sottoarea_splittata[5]))
    if tot_area==0:
        avrg_density=0
    else:
        avrg_density=tot_objects/tot_area
    if stdevs_list==[]:
        avrg_stdev=0
    else:
        avrg_stdev=stat.mean(stdevs_list)
    area_unita=sigla+';'+nome_completo+';'+str(tot_objects)+';'+str(tot_area)+';'+str(avrg_density)+';'+str(avrg_stdev)+';'+'0'+';'+'0'+'\n'
    return area_unita
    
    
def unite_special(lista_righe, first_line, last_line, additional_line, sigla, nome_completo):
    tot_objects=0
    tot_area=0
    avrg_density=0
    stdevs_list=[]
    avrg_stdev=0
    for i in range(first_line, last_line+1):
        sottoarea_splittata=lista_righe[i].split(';')
        tot_objects+=float(sottoarea_splittata[2])
        tot_area+=float(sottoarea_splittata[3])
        stdevs_list.append(float(sottoarea_splittata[5]))
    sottoarea_splittata=lista_righe[additional_line].split(';')
    tot_objects+=float(sottoarea_splittata[2])
    tot_area+=float(sottoarea_splittata[3])
    stdevs_list.append(float(sottoarea_splittata[5]))
    if tot_area==0:
        avrg_density=0
    else:
        avrg_density=tot_objects/tot_area
    if stdevs_list==[]:
        avrg_stdev=0
    else:
        avrg_stdev=stat.mean(stdevs_list)
    area_unita=sigla+';'+nome_completo+';'+str(tot_objects)+';'+str(tot_area)+';'+str(avrg_density)+';'+str(avrg_stdev)+';'+'0'+';'+'0'+'\n'
    return area_unita


def copy(lista_righe, first_line, last_line):
    lista_righe_da_copiare=lista_righe[first_line:last_line+1]
    return lista_righe_da_copiare

with open('aree_da_unire_eliminare_copiare_per_scendere_a_134.txt', 'r') as f:
    lista_istruzioni=f.readlines()
lista_file=get_list_of_file_in_directory('con_sigle')

for file in lista_file:
    with open(file, 'r') as r, open('134_strutture/134_'+file.split('\\')[1], 'a') as w:
        lista_righe=r.readlines()
        for i in lista_istruzioni:
            istruzione=i.split(',')
            first_line=int(istruzione[1])
            last_line=int(istruzione[2].rstrip('\n'))
            if istruzione[0]=='E':
                pass
            elif istruzione[0]=='C':
                w.writelines(copy(lista_righe, first_line, last_line))
            elif istruzione[0]=='U':
                print(istruzione)
                w.write(unite(lista_righe, first_line, last_line, istruzione[3], istruzione[4].rstrip('\n')))
            elif istruzione[0]=='US':
                w.write(unite_special(lista_righe, first_line, last_line, int(istruzione[3]), istruzione[4], istruzione[5].rstrip('\n')))
            elif istruzione[0]=='R':
                w.write(rename(lista_righe, first_line, last_line, istruzione[3], istruzione[4].rstrip('\n')))
            else:
                print(istruzione)
                print("ERRORE")



            
    
    
    
    
    
    
