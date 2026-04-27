import os
import statistics

def get_list_of_files_in_directory(dir): #ottiene lista file in una directory escludendo le altre directory ivi contenute (e in generale tutto ciò che non è file)
    if(not os.path.isdir(dir)):
        raise Exception(f"{dir} is not a directory")
    files = []
    for f in os.scandir(dir):
        if(f.is_file()):
            files.append(f.path)
    return files

#tocca fa un dizionario con chiave = numero fettina e contenuto = lista annidata di liste (a.k.a. righe splittate)
def create_dict_Key_csv_Value_lista_annidata_righe (lista_csv):
    dictionary={}
    for csv in lista_csv:
        if csv[-11:]=='Regions.csv':
            pass
        else:
            lista_annidata_righe=[]
            with open(csv, 'r') as f:
                lista_annidata_righe = f.readlines()
            counter = 0
            for i in lista_annidata_righe:
                lista_annidata_righe[counter]=lista_annidata_righe[counter].split(';')
                counter += 1
            dictionary[csv]=lista_annidata_righe
    return dictionary

#STA ROBA COMMENTATA NON SERVE PERCHé SONO SOLO DUE ID E DUE COUNT E GLIELI DO IO
# def obtain_lista_id_e_dizionario(lista_id, lista_count_ideali):
#     dizionario={}
#     with open(lista_id, 'r') as f:
#         lista_id_provvisoria=f.read()
#     lista_id_definitiva=lista_id_provvisoria.split('\n')
#     with open(lista_count_ideali, 'r') as f:
#         lista_count_ideali_provvisoria=f.read()
#     lista_count_ideali_definitiva=lista_count_ideali_provvisoria.split('\n')
#     if len(lista_count_ideali_definitiva) == len(lista_id_definitiva):
#         for i in range(0, len(lista_id_definitiva)):
#             dizionario[lista_id_definitiva[i]]=lista_count_ideali_definitiva[i]
#     else:
#         with open('ERRORE.txt', 'w') as f:
#             f.write('AAAAAAAAAA!!!!!!!!!!!!!!!!!!')
#     return lista_id_definitiva, dizionario



def obtain_id_name_totcount_totarea_avrgdensity_countfettine_countideale(dictionary, lista_id, dizionario_idVScount):
    listarighe_nuovofile=[]
    for id in lista_id:
        name=''
        totcount=0
        totarea=0
        densities_list=[]
        avrgdensity=0
        densities_stdev=0
        countfettine=0
        countideale=dizionario_idVScount[id]
        for key in dictionary.keys():
            for riga_splittata in dictionary[key]:
                if riga_splittata[0]==id and id=='672':
                    if riga_splittata[2]=='0' or riga_splittata[2]=='NA':
                        name = 'DLS'
                    else:
                        name = 'DLS'
                        totcount += float(riga_splittata[5])
                        totarea += float(riga_splittata[11])
                        densities_list.append(float(riga_splittata[12]))
                        countfettine += 1
            for riga_splittata in dictionary[key]:
                if riga_splittata[0]==id and id=='56':
                    if riga_splittata[2]=='0' or riga_splittata[2]=='NA':
                        name = 'NACC_SHELL'
                    else:
                        name = 'NACC_SHELL'
                        totcount += float(riga_splittata[5])
                        totarea += float(riga_splittata[11])
                        densities_list.append(float(riga_splittata[12]))
                        countfettine += 1
        if countfettine>0:
            print(densities_list)
            avrgdensity=totcount/totarea
            densities_stdev=statistics.pstdev(densities_list)
        listarighe_nuovofile.append(id+'tris;'+name+';'+str(totcount)+';'+str(totarea)+';'+str(avrgdensity)+';'+str(densities_stdev)+';'+str(countfettine)+';'+str(countideale)+'\n')
    return listarighe_nuovofile


lista_csv = get_list_of_files_in_directory('MASK_RefAtlas_with_density')
dictionary = create_dict_Key_csv_Value_lista_annidata_righe(lista_csv)
#lista_id, dizionario_idVScount = obtain_lista_id_e_dizionario('id_da_mettere.txt', 'count_ideali.txt')
#caudoputamen ha ID 672 e count ideale 39, nucleus accumbens ha ID 56 e count ideale 16
lista_id=['672', '56']
dizionario_idVScount={}
dizionario_idVScount['672']=39
dizionario_idVScount['56']=16
with open('MASK_'+os.getcwd()[-3:]+'_unito.csv', 'a')as f:
    f.write('id;name;tot_objects;tot_area;avrg_density;densities_stdev;slices_count;ideal_count\n')
    f.writelines(obtain_id_name_totcount_totarea_avrgdensity_countfettine_countideale(dictionary, lista_id, dizionario_idVScount))

