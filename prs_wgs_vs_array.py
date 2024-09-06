# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 19:35:01 2024

@author: Edurne Urrutia-Lafuente

 PRS: WGS vs array comparison © 03.2024 by Edurne Urrutia-Lafuente is licensed
 under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.
 To view a copy of this license, visit
 https://creativecommons.org/licenses/by-nc-sa/4.0/
"""

import pandas as pd
import csv

# Paso 1: Cargar el archivo VCF multisample
# Definir una lista para almacenar los datos
vcf_data = []
smpl_lst = smpl_lst.split("', '")
snp_file = "C:/Users/D822886/Desktop/MX/CHIP/chip_snps.csv"

# Inicializar una lista para almacenar los datos del archivo CSV
snp_lst = []

# # Abrir el archivo CSV y leer los datos en la lista
# with open(snp_file, newline='') as csvfile:
#     lector_csv = csv.reader(csvfile)
#     for fila in lector_csv:
#         snp_lst.append(fila)


# Abrir el archivo CSV y leer los datos en la lista
with open(snp_file, newline='') as csvfile:
    lector_csv = csv.reader(csvfile)
    for fila in lector_csv:
        # Convertir cada elemento de la fila a string y unirlos en un solo string
        fila_str = ','.join(fila)
        snp_lst.append(fila_str)


# Abrir el archivo y leerlo línea por línea
with open(C:/Users/D822886/Desktop/MX/CHIP/chip_snps.csvMuestras.hg19.313 (2).vcf", "r") as f:
    for line in f:
        # Ignorar las líneas que comienzan con ##
        if line.startswith("##"):
            continue
        if line.startswith('#'):
            # Si es una línea que comienza con # pero no con ##, es el header, podemos obtener los nombres de las columnas
            header = line.strip().lstrip("#").split("\t")

        # Procesar solo las líneas que no empiezan con '#' (comentarios)
        elif not line.startswith("#"):
            # Separar la línea por tabulaciones y agregar los datos a la lista
            row = line.strip().split("\t")
            vcf_data.append(row)


# Convertir la lista en un DataFrame de pandas
vcf_df = pd.DataFrame(vcf_data, columns=header)

# Quitar el prefijo "chr" de la columna "Chromosome"
vcf_df["CHROM"] = vcf_df["CHROM"].str.replace("chr", "")

# Combinar las columnas "Chromosome" y "Position" para formar el índice en vcf_df
vcf_df.index = vcf_df["CHROM"].astype(str) + "_" + vcf_df["POS"].astype(str)

#smpl list completa del vcf
column_names = vcf_df.columns.tolist()[9:]

for name in column_names:
    if name not in smpl_lst:
        print(name)
        
        
        
# Paso 2: Inicializar diccionario para almacenar las coincidencias por posición
coincidencias_por_muestra = {}
coincidencias_por_posicion = {}
# Inicializar un diccionario para almacenar las posiciones faltantes y el número de muestras en las que no están
posiciones_faltantes = {}

# Iniciar diccionarios
glob_dct = {"total":0, "ok":0, "mal":0, "nan":0}
pos_dct = {}
smpl_dct = {}
msimatches = []

# Paso 3: Iterar sobre cada archivo TXT de secuenciación
# Cargar el archivo de equivalencias
equivalencias_df = pd.read_csv("C:/Users/D822886/Desktop/ref_cnag.csv", sep = ';')

# Inicializar un diccionario para mapear referencias de CNAG a referencias de Mx
cnag_a_mx = dict(zip(equivalencias_df["RefNageMX"], equivalencias_df["ReferenciaCNAG"]))



# Identificar filas duplicadas
filas_duplicadas = vcf_df.index.duplicated(keep='first')

# Mantener solo la primera ocurrencia de las filas duplicadas
vcf_df_uniq = vcf_df[~filas_duplicadas]

for smplmx in smpl_lst:
    smpl = cnag_a_mx[smplmx]
    txt_file = f'T:/NAGEN-Mx/ANÁLISIS/1_PRS_NAGEN-Mx/0_SCRIPT_FILES/Todos/{smpl}.PRS_313_SNPs.DPge10.txt'
    
    wgs_df = pd.read_csv(txt_file, sep="\t")
    
    # Combinar las columnas "Chromosome" y "Position" para formar el índice
    wgs_df.index = wgs_df["Chromosome_x"].astype(str) + "_" + wgs_df["Position"].astype(str)
    
    # Iterar sobre las posiciones de la lista:
    for snp in snp_lst:
        if snp in wgs_df.index:
            wgs_gt = wgs_df.at[snp, "NumGT"]
            if snp == '10_123340431':
                chip_gt = str(vcf_df_uniq.at['10_123340432',smplmx]).split(":")[0]
            else:
                chip_gt = str(vcf_df_uniq.at[snp,smplmx]).split(":")[0]
            if wgs_gt == chip_gt:
                glob_dct["ok"] += 1
            else:
                if chip_gt == './.':
                    glob_dct["nan"] += 1
                    mismatches.append([smplmx, smpl, snp, wgs_gt, chip_gt])
                elif snp == '7_99948655':
                    wgs_gt_new = wgs_gt.replace("2", "0")
                    if wgs_gt_new == chip_gt:
                        glob_dct["ok"] += 1
                    else:
                        print(smplmx, smpl, snp, wgs_gt_new, chip_gt)
                        glob_dct["mal"] += 1
                        mismatches.append([smplmx, smpl, snp, wgs_gt_new, chip_gt])
                elif snp == '22_40904707':
                    wgs_gt_new = wgs_gt.replace("2", "0")
                    if wgs_gt_new == chip_gt:
                        glob_dct["ok"] += 1
                    else:
                        print(smplmx, smpl, snp, wgs_gt_new, chip_gt)
                        glob_dct["mal"] += 1
                        mismatches.append([smplmx, smpl, snp, wgs_gt_new, chip_gt])
                elif snp == '3_141112859':
                    wgs_gt_new = wgs_gt.replace("2", "0")
                    if wgs_gt_new == chip_gt:
                        glob_dct["ok"] += 1
                    else:
                        print(smplmx, smpl, snp, wgs_gt_new, chip_gt)
                        glob_dct["mal"] += 1
                        mismatches.append([smplmx, smpl, snp, wgs_gt_new, chip_gt])
                elif snp == '5_52679539':
                    wgs_gt_new = wgs_gt.replace("2", "0")
                    if wgs_gt_new == chip_gt:
                        glob_dct["ok"] += 1
                    else:
                        print(smplmx, smpl, snp, wgs_gt_new, chip_gt)
                        glob_dct["mal"] += 1
                        mismatches.append([smplmx, smpl, snp, wgs_gt_new, chip_gt])
                elif snp == '15_100905819':
                    wgs_gt_new = wgs_gt.replace("2", "0")
                    if wgs_gt_new == chip_gt:
                        glob_dct["ok"] += 1
                    else:
                        print(smplmx, smpl, snp, wgs_gt_new, chip_gt)
                        glob_dct["mal"] += 1
                        mismatches.append([smplmx, smpl, snp, wgs_gt_new, chip_gt])
                else:
                    print(smplmx, smpl, snp, wgs_gt, chip_gt)
                    glob_dct["mal"] += 1
                    mismatches.append([smplmx, smpl, snp, wgs_gt, chip_gt])
        else:
            glob_dct["nan"] += 1
            #pos_dct

# Crear un DataFrame a partir de la lista de no coincidencias
mismatch_df = pd.DataFrame(mismatches, columns=["Ref. NagenMX", "Ref. CNAG", "Posicion", "Genotipo WGS", "Genotipo CHIP"])

# Escribir el DataFrame en un archivo Excel
mismatch_df.to_excel("C:/Users/D822886/Desktop/MX/CHIP/mismatches.xlsx", index=False)
        
    #for index, row in wgs_df.iterrows():
        chrom = row["Chromosome_x"]
        position = row["Position"]
        locus = index

        # Verificar si la posición está presente en el archivo VCF
        if index not in vcf_df.index:
            posiciones_faltantes.setdefault(index, 0)
            posiciones_faltantes[index] += 1
            continue
        
        # Obtener el genotipo correspondiente desde el archivo VCF
        vcf_genotipo = str(vcf_df.loc[index, smplmx]).split(':')[0]
        #print(vcf_genotipo)
        #print(index, smplmx)
        #print(vcf_genotipo)
        #print(vcf_df.loc[index, smplmx])
        # Obtener el genotipo del archivo TXT
        print(row["NumGT"])
        wgs_genotipo = row["NumGT"].split(':')[0]
    
        # Paso 4: Comparar los genotipos del archivo TXT con los del archivo VCF
        if vcf_genotipo
        coincidencias = (vcf_df["GT"] == wgs_df["GT"]).astype(int)
    
    # Paso 5: Actualizar el diccionario de coincidencias por posición
    for i, position in enumerate(vcf_df["Position"]):
        coincidencias_por_posicion.setdefault(position, {})
        coincidencias_por_posicion[position][smpl] = coincidencias[i]

# Paso 6: Calcular el número total de coincidencias globalmente
coincidencia_global = sum(sum(sample.values()) for sample in coincidencias_por_posicion.values())

# Mostrar los resultados
print("Coincidencia global:", coincidencia_global)
print("Coincidencia por posición:")
for position, samples in coincidencias_por_posicion.items():
    print(f"Posición {position}: {sum(samples.values())} coincidencias")
    print("   Muestras:", samples)