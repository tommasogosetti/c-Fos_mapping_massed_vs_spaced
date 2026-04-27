import statistics as stat
import os

#COSE DA FARE - eliminare riga 2, fare righe nel file di testo, eliminare tutto ciò che c'è da 529, unire 320-324 e 325-329

#in questa funzione il count fettine dell'area unita è la media delle singole sottoaree mentre il count ideale, non avendo senso, è 0, ricordarsi di cambiarlo se non andasse bene
#la variabile avrg_countfettine è calcolata ma non inserita, perché il prof ha detto di non metterlo: per integrarlo basta sostituire lo 0 con la variabile
def unisci_sottoaree(lista_righe, area, prima_sottoarea, ultima_sottoarea):
    area_splittata=lista_righe[area-1].split(';')
    tot_objects=0
    tot_area=0
    densities_list=[]
    avrg_density=0
    stdevs_list=[]
    avrg_stdev=0
    countfettine_list=[]
    avrg_countfettine=0
    countideale=0
    for i in range(prima_sottoarea-1, ultima_sottoarea):
        sottoarea_splittata=lista_righe[i].split(';')
        if sottoarea_splittata[3]!='0':
            tot_objects+=float(sottoarea_splittata[2])
            tot_area+=float(sottoarea_splittata[3])
            densities_list.append(float(sottoarea_splittata[4]))
            stdevs_list.append(float(sottoarea_splittata[5]))
            countfettine_list.append(float(sottoarea_splittata[6]))
    if tot_area==0:
        avrg_density=0
    else:
        avrg_density=tot_objects/tot_area
    if stdevs_list==[]:
        avrg_stdev=0
    else:
        avrg_stdev=stat.mean(stdevs_list)
    if countfettine_list==[]:
        avrg_countfettine=0
    else:
        avrg_countfettine=stat.mean(countfettine_list)
#    avrg_density, avrg_stdev, avrg_countfettine = stat.mean(densities_list), stat.mean(stdevs_list), stat.mean(countfettine_list)
    area_unita=area_splittata[0]+';'+area_splittata[1]+';'+str(tot_objects)+';'+str(tot_area)+';'+str(avrg_density)+';'+str(avrg_stdev)+';'+'0'+';'+str(countideale)+'\n'
    lista_righe[area-1]=area_unita
    lista_righe[prima_sottoarea-1:ultima_sottoarea]=[]
    return lista_righe


# da qui in poi ESEGUIRE FUNZIONI NELL'ORDINE IN CUI LE DEFINISCO
#(sennò si hanno problemi nella numerazione delle righe)


#elimina varie fibre
def elimina_righe_da_529esima(csv):
    with open(csv, 'r') as f:
        lista_righe=f.readlines()
    lista_righe[528:]=[]
    return lista_righe

# a livello teoricose le sottoaree sono solo 2 potrebbe dare problemi,
# ma non è il nostro caso dato che dobbiamo unire 320-324 e 325-329.
# idem nel caso in cui le righe fossero sia prima che dopo quelle manipolate da altre funzioni
#(poiché cambieremmo gli indici della lista righe), ma non è il nostro caso.
#se un giorno volessi toccare anche il countideale, ricordarsi che c'è un \n alla fine di ogni riga e va rstrippato
#anche qui la variabile avrg_countfettine è calcolata ma non inserita, perché il prof ha detto di non metterlo: per integrarlo basta sostituire lo 0 con la variabile
def unisci_sottoaree_senza_area(lista_righe, id_area, nome_area, prima_sottoarea, ultima_sottoarea):
    tot_objects=0
    tot_area=0
    densities_list=[]
    avrg_density=0
    stdevs_list=[]
    avrg_stdev=0
    countfettine_list=[]
    avrg_countfettine=0
    countideale=0
    for i in range(prima_sottoarea-1, ultima_sottoarea):
        sottoarea_splittata=lista_righe[i].split(';')
        if sottoarea_splittata[3]!='0':
            tot_objects+=float(sottoarea_splittata[2])
            tot_area+=float(sottoarea_splittata[3])
            densities_list.append(float(sottoarea_splittata[4]))
            stdevs_list.append(float(sottoarea_splittata[5]))
            countfettine_list.append(float(sottoarea_splittata[6]))
    if tot_area==0:
        avrg_density=0
    else:
        avrg_density=tot_objects/tot_area
    if stdevs_list==[]:
        avrg_stdev=0
    else:
        avrg_stdev=stat.mean(stdevs_list)
    if countfettine_list==[]:
        avrg_countfettine=0
    else:
        avrg_countfettine=stat.mean(countfettine_list)
#    avrg_density, avrg_stdev, avrg_countfettine = stat.mean(densities_list), stat.mean(stdevs_list), stat.mean(countfettine_list)
    area_unita=id_area+';'+nome_area+';'+str(tot_objects)+';'+str(tot_area)+';'+str(avrg_density)+';'+str(avrg_stdev)+';'+'0'+';'+str(countideale)+'\n'
    lista_righe[prima_sottoarea-1]=area_unita
    lista_righe[prima_sottoarea:ultima_sottoarea]=[]
    return lista_righe

#ciclo for reversato in maniera da andare dalle ultime righe alle prime e non confondere l'ordine
def unisci_sottoaree_da_lista(lista_righe, sottolayer_da_scremare):
    with open(sottolayer_da_scremare, 'r') as f:
        lista_sottolayer=f.readlines()
        print(lista_sottolayer)
    for i in reversed(lista_sottolayer):
        j=i.rstrip('\n')
        j=j.split(',')
        area=int(j[0])
        prima_sottoarea=int(j[1])
        ultima_sottoarea=int(j[2])
        lista_righe = unisci_sottoaree(lista_righe, area, prima_sottoarea, ultima_sottoarea)
    return lista_righe


#elimina root
def elimina_riga_2(lista_righe):
    lista_righe.pop(1)
    return lista_righe


#elimina orbital ventral area, lasciando solo orbital ventrolateral area (con valori medie e somme di quelli dei suoi sottolayer)
def elimina_riga_33(lista_righe):
    lista_righe.pop(32)
    return lista_righe


lista_righe = elimina_righe_da_529esima(os.getcwd()[-3:]+'_unito.csv')
lista_righe = unisci_sottoaree_senza_area(lista_righe, '526,543,664,727,743', 'Entorhinal area, medial part, dorsal zone', 325, 329)
lista_righe = unisci_sottoaree_senza_area(lista_righe, '1121,20,52,139,28', 'Entorhinal area, lateral part', 320, 324)
lista_righe = unisci_sottoaree_da_lista(lista_righe, 'sottolayer_da_scremare.txt')
lista_righe = elimina_riga_2(lista_righe)
lista_righe = elimina_riga_33(lista_righe)
with open(os.getcwd()[-3:]+'_scremato.csv', 'w') as f:
    f.writelines(lista_righe)
