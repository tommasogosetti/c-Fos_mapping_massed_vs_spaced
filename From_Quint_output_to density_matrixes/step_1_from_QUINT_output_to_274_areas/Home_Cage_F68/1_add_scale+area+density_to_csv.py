import os

def get_list_of_file_in_directory(dir): #ottiene lista file in una directory escludendo le altre directory ivi contenute (e in generale tutto ciò che non è file)
    if(not os.path.isdir(dir)):
        raise Exception(f"{dir} is not a directory")
    files = []
    for f in os.scandir(dir):
        if(f.is_file()):
            files.append(f.path)
    return files

def create_lista_annidata_fettine (results):
    lista_annidata_fettine=[]
    counter = 0
    with open(results, 'r') as f:
        lista_immagini=f.readlines()
    for fettina in lista_immagini:
        if counter == 0:
            counter = 1
        elif counter == 1:
            sublist=fettina.split(', ')
            lista_annidata_fettine.append(sublist)
    return lista_annidata_fettine


def create_dict (lista_annidata_fettine, lista_file_csv):
    dict_fettinaVScsv={}
    for fetta in lista_annidata_fettine:
        for csv in lista_file_csv:
            if(fetta[0])[-7:-4]==csv[-7:-4]:
                dict_fettinaVScsv[csv]=fetta[1]
    return dict_fettinaVScsv

# def add_scale_and_density_to_every_csv _VECCHIO (dict_fettinaVScsv, lista_file_csv):
#     for csv in lista_file_csv:
#         with open(csv, 'r') as f:
#             numero_righe=len(f.readlines())
#         with open(csv, 'r+') as f:
#             for l in numero_righe:
#                 riga=f.readline()
#                 riga_splittata=riga.split(';')
#                 if int(riga_splittata[2])!=0:
#                     area=float(riga_splittata[2])*float(dict_fettinaVScsv[csv])*float(dict_fettinaVScsv[csv])
#                     density=float(riga_splittata[5])/area
#                     f.write(';'+dict_fettinaVScsv[csv]+';'+str(area)+';'+str(density)+'\n')
                    
def add_scale_and_density_to_every_csv (dict_fettinaVScsv, cartella):
    for csv in list(dict_fettinaVScsv.keys()):
        with open(csv, 'r') as f:
            listone_righe = f.readlines()
            counter = 0
        for riga in listone_righe:
            #questo toglie le colonne, cioè i separatori, in eccesso dalla variabile "riga"
            riga=riga.rstrip('\n')
            while riga[-1]==';':
                riga=riga[:-1]
            riga+='\n'
            #fine parte che toglie i separatori in eccesso, inizio parte che, solo per la prima riga, corregge l'header
            if counter == 0:
                riga_corretta=riga.rstrip('\n')
                riga_corretta+=(';Scale(micron/pixel);Area(micron2);Density(cells/micron2)\n')
                listone_righe[counter]=riga_corretta
                print(listone_righe[0])
                counter = 1
            #fine parte dedicata all'header, inizio parte in cui effettivamente aggiunge i dati alle righe successive
            elif counter > 0:
                riga_splittata = riga.split(';')
                try:
                    if int(riga_splittata[2])!=0:
                        area=float(riga_splittata[2])*float(dict_fettinaVScsv[csv])*float(dict_fettinaVScsv[csv])
                        density=float(riga_splittata[5])/area
                        riga_corretta=riga.rstrip('\n')
                        riga_corretta+=(';'+dict_fettinaVScsv[csv]+';'+str(area)+';'+str(density)+'\n')
                        listone_righe[counter]=riga_corretta
                        print(listone_righe[counter])
                        counter+=1
                    else:
                        listone_righe[counter]=riga
                        counter+=1
                #except necessario per gli N/A che non possono essere letti come int
                except ValueError:
                    listone_righe[counter]=riga
                    counter+=1
        with open(cartella+'/'+os.getcwd()[-3:]+csv[-10:], 'a') as f:
            f.writelines(listone_righe)
    
                   

lista_file_csv = get_list_of_file_in_directory('RefAtlasRegions')
lista_annidata_fettine = create_lista_annidata_fettine('scales.txt')
print(lista_annidata_fettine)
dict_fettinaVScsv = create_dict(lista_annidata_fettine, lista_file_csv)
cartella='RefAtlas_with_density'
os.mkdir(cartella)
add_scale_and_density_to_every_csv (dict_fettinaVScsv, cartella)
print(os.getcwd())