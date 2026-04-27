import os

def ottieni_listarighe_normale_e_mask (csv_normale, csv_mask):
    with open(csv_normale, 'r') as normale, open(csv_mask, 'r') as mask:
        listarighe_normale=normale.readlines()
        listarighe_mask=mask.readlines()
    return listarighe_normale, listarighe_mask

def ottieni_DLS_DMS_NaccShell_NaccCore (listarighe_normale, listarighe_mask):
    dls=listarighe_mask[1]
    nac_shell=listarighe_mask[2]
    caudoputamen_spl=listarighe_normale[89].split(';')
    dls_spl=dls.split(';')
    dms='672bis;DMS;'+str(float(caudoputamen_spl[2])-float(dls_spl[2]))+';'+str(float(caudoputamen_spl[3])-float(dls_spl[3]))+';'+str((float(caudoputamen_spl[2])-float(dls_spl[2]))/(float(caudoputamen_spl[3])-float(dls_spl[3])))+';;'+caudoputamen_spl[6]+';'+caudoputamen_spl[7]
    nac_spl=listarighe_normale[90].split(';')
    nac_shell_spl=nac_shell.split(';')
    nac_core='56bis;NACC_CORE;'+str(float(nac_spl[2])-float(nac_shell_spl[2]))+';'+str(float(nac_spl[3])-float(nac_shell_spl[3]))+';'+str((float(nac_spl[2])-float(nac_shell_spl[2]))/(float(nac_spl[3])-float(nac_shell_spl[3])))+';;'+nac_spl[6]+';'+nac_spl[7]
    return dls, dms, nac_shell, nac_core

def inserisci_valori_nella_lista (dls, dms, nac_shell, nac_core, listarighe_normale):
    listarighe_normale.insert(90, nac_core)
    listarighe_normale.insert(90, nac_shell)
    listarighe_normale.insert(89, dms)
    listarighe_normale.insert(89, dls)
    return listarighe_normale



listarighe_normale, listarighe_mask = ottieni_listarighe_normale_e_mask(os.getcwd()[-3:]+'_scremato.csv', 'MASK_'+os.getcwd()[-3:]+'_unito.csv')
dls, dms, nac_shell, nac_core = ottieni_DLS_DMS_NaccShell_NaccCore(listarighe_normale, listarighe_mask)
listarighe_finale = inserisci_valori_nella_lista(dls, dms, nac_shell, nac_core, listarighe_normale)
with open(os.getcwd()[-3:]+'_final_for_variability.csv', 'w') as f:
    f.writelines(listarighe_finale)
