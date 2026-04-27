import os

def lista_topi_gruppo(sigla_gruppo):
    gruppo=[]
    for i in os.listdir(os.getcwd()):
        print(i)
        if sigla_gruppo in i:
            gruppo.append(i)
    return gruppo

            
HC=lista_topi_gruppo('_HC_')
MS=lista_topi_gruppo('_MS_')
MC=lista_topi_gruppo('_MC_')
SS=lista_topi_gruppo('_SS_')
SC=lista_topi_gruppo('_SC_')


def obtain_header():
    with open('134_SIGLE_HC_F39_final_for_variability.csv', 'r') as f: #f39 è uno qualunque. tanto l'header è uguale per tutti
        listarighe=f.readlines()
    header='animale;'
    prima_riga='zompa'
    for i in listarighe:
        if prima_riga=='zompa':
            prima_riga='passata'
        elif prima_riga=='passata':
            riga_splittata=i.split(';')
            if riga_splittata[0] in lista_aree_da_eliminare:
                print(riga_splittata[0]+' trovata e non inserita in header')
            else:
                header+=riga_splittata[0]+';'
        else:
            print("ERRORE")
    return header

def obtain_animal(animale, lista_aree_da_eliminare, lista_aree_da_tenere):
    with open(animale, 'r') as f:
        listarighe=f.readlines()
    nuova_riga=animale[13:-26]+';'
    prima_riga='zompa'
    for i in listarighe:
        if prima_riga=='zompa':
            prima_riga='passata'
        elif prima_riga=='passata':
            riga_splittata=i.split(';')
            if riga_splittata[0] in lista_aree_da_eliminare:
                print(riga_splittata[0]+' trovata e non inserita!')
            elif riga_splittata[6]=='1' or riga_splittata[6]=='2' and riga_splittata[0] not in lista_aree_da_tenere:
                nuova_riga+='NA;'
            elif riga_splittata[6]=='0' and riga_splittata[4]=='0':
                nuova_riga+='NA;'
            else:
                nuova_riga+=riga_splittata[4]+';'
        else:
            print("ERRORE")
    return nuova_riga

def generate_group_file_content(group, header):
    content=[header]
    for i in group:
        riga_animale=obtain_animal(i, lista_aree_da_eliminare, lista_aree_da_tenere)+'\n'
        content.append(riga_animale)
    return content


lista_aree_da_eliminare=['parVPT', 'SGn', 'ADn', 'LDT', 'VlPO', 'LMN', 'lMMN', 'VPrM'] #MODIFICABILE SE DOVESSIMO CAMBIARE AREE
lista_aree_da_tenere=['cSN', 'VTA', 'reSn', 'Gpi', 'Gpe']

header=obtain_header()+'\n'

with open('file_dei_gruppi\HC_per_matrici.csv', 'w') as f:
    f.writelines(generate_group_file_content(HC, header))
    
with open('file_dei_gruppi\MC_per_matrici.csv', 'w') as f:
    f.writelines(generate_group_file_content(MC, header))
    
with open('file_dei_gruppi\SC_per_matrici.csv', 'w') as f:
    f.writelines(generate_group_file_content(SC, header))
    
with open('file_dei_gruppi\SS_per_matrici.csv', 'w') as f:
    f.writelines(generate_group_file_content(SS, header))
    
with open('file_dei_gruppi\MS_per_matrici.csv', 'w') as f:
    f.writelines(generate_group_file_content(MS, header))