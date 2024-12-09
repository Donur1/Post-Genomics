# import pandas as pd
# import pandas as pd
# from lifelines import CoxPHFitter

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from lifelines import CoxPHFitter
#


#
# ############ VCF to CSV ############
# def splitINFO(file): # file would have to be a dataframe
#     csq = []
#     for row in range(file.shape[0]):
#         info1 = str(file["info"][row])
#         start = info1.find("CSQ")
#         end = info1.find("SOMATIC")
#         if end > start:
#             info = info1[start+6:end-3]
#             info = info.replace("'", "")
#             info = info.replace("[", "")
#             info = info.replace("]", "")
#             info = info.replace(" ", "")
#             info = info.split(",")
#         else:
#             info = info1[start+6:-2]
#             info = info.replace("'", "")
#             info = info.replace("[", "")
#             info = info.replace("]", "")
#             info = info.replace(" ", "")
#             info = info.split(",")
#
#         csq.append(info)
#     file.insert(file.shape[1], "CSQ", csq)
#     return file
#
#
#
# ################### separate CSQ list into individual elements for new columns ######################
#     # AML headers
# AMLheader = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|ENTREZ|EVIDENCE"
#     # ALL headers
# ALLheader = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|ENTREZ|EVIDENCE"
#     # Thyroid headers
# THCAheader = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS"
#     # Ovarian headers
# thcaheader = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS"
#
#
# def parseCSQ (file, header): #file needs to be a dataframe, header needs to be based on the information within the vcf
#     newColheader = header.split("|")
#
#     # creates new column headers based on inputted list where the entries are empty list that can be appended to:
#     for col in newColheader:
#         file[col] = [list() for x in range(len(file.index))]
#
#
#     for row in range(file.shape[0]):
#         info1 = str(file["CSQ"][row])
#         info1 = info1.split(",")
#
#         for ele in info1:
#             count = 0
#             ele = ele.split("|")
#
#             while count < len(ele):
#                 file[newColheader[count]][row].append(ele[count])
#                 count += 1
#     return file
#
# def CaseIDs(vcffile):
#     temp = vcffile.split(".")
#     vcffile = open('%s'%(vcffile), 'r')
#     vcf = vcffile.readlines()
#     vcffile.close()
#
#     for l in range(0, len(vcf)): #looking line by line in the inputted VCF file
#         if vcf[l].startswith('##INDIVIDUAL'):
#
#             ind1 = vcf[l].find("ID=")
#             ind2 = vcf[l].find("NAME=")
#
#             line = vcf[l]
#             if ind1 == -1:
#                 print("Problem with: " + temp[0])
#                 return
#             else:
#                 #print(line[ind+3:-2])
#                 ele1 = temp[1] #Patient ID/VCF file name
#                 ele2 = line[ind1+3:-2] #Case ID that will be matched with cBio portal
#                 ele3 = line[ind2+5:ind1-1]
#                 #print(row)
#
#     return ele1, ele2, ele3
#
# def addALT_Seq(csv):
#     alt = []
#     for row in range(csv.shape[0]):
#         ref_seq = csv["ref_seq"][row]
#         if ref_seq == csv["var_seq1"][row]:
#             alt.append(csv["var_seq2"][row])
#         else:
#             alt.append(csv["var_seq1"][row])
#     csv.insert(csv.shape[1], "alt_seq", alt)
#     return csv
#
# def writecsv(inputfile, header):
#     vcf_reader = vcf.Reader(open(inputfile,'r'))  # read vcf file
#     #use nested list to store variant information, one for tumor, one for normal
#     normal_var=[]
#     tumor_var=[]
#     var_score = [28]
#     refct = []
#     altct = []
#     refct1 = []
#     altct1 = []
#
#     for record in vcf_reader:
#         curr_var = []
#         left = []
#         right = []
#         ref_seq = []
#         alt_seq = []
#         info = []
#         # Because the AML and ALL have a reverse order for which comes first in the VCF, the sample[0/1] needed to be changed in order to make sure
#             # the tumor and normal variants were saved to the correct csv. Can be compared to the original PrCa code
#         tumorAD = record.samples[1]["AD"]
#         normalAD = record.samples[0]["AD"]
#         ########################################### tumor_var #############################################
#         #######check for 2 ALT options:
#         if len(tumorAD) == 3:
#             sum = tumorAD[0] + tumorAD[1] + tumorAD[2]
#             x = tumorAD[0]
#             y = tumorAD[1]
#             z = tumorAD[2]
#             # if two out of the three are not greater than 0.05 throw away variant
#             if ((x / sum) < 0.05 and (y / sum) < 0.05) or ((x / sum) < 0.05 and (z / sum) < 0.05) or (
#                     (z / sum) < 0.05 and (y / sum) < 0.05):
#                 pass
#             # means two are going to be good to decide which to use
#             else:
#                 print(record)
#
#         ########only 1 base for the ALT sequence
#         else:
#             sum = tumorAD[0] + tumorAD[1]
#             refcount = tumorAD[0]
#             altcount = tumorAD[1]
#
#             # both the REF and ALT have to be greater
#             if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
#                 curr_var.append(record.CHROM)
#                 left.append(record.POS)
#                 right.append(record.POS + len(record.ALT))
#                 ref_seq.append(record.REF)
#                 alt_seq.append(str((record.ALT[0])))
#                 refct.append(refcount)
#                 altct.append(altcount)
#                 info.append(str(record.INFO))
#                 lists = curr_var + left + right + ref_seq + alt_seq + var_score + info
#                 tumor_var.append(lists)
#
#
#             # throw away variant if not
#             else:
#                 pass
#
#         ############################# normal_var ##############################3
#         curr_var1 = []
#         left1 = []
#         right1 = []
#         ref_seq1 = []
#         alt_seq1 = []
#         info = []
#         if len(normalAD) == 3:
#             sum = normalAD[0] + normalAD[1] + normalAD[2]
#             x = normalAD[0]
#             y = normalAD[1]
#             z = normalAD[2]
#             # if two out of the three are not greater than 0.05 throw away variant
#             if ((x / sum) < 0.05 and (y / sum) < 0.05) or ((x / sum) < 0.05 and (z / sum) < 0.05) or (
#                     (z / sum) < 0.05 and (y / sum) < 0.05):
#                 pass
#             # means two are going to be good to decide which to use
#             else:
#                 print(record)
#
#         ########only 1 base for the ALT sequence
#         else:
#             sum = normalAD[0] + normalAD[1]
#             refcount = normalAD[0]
#             altcount = normalAD[1]
#
#             # both the REF and ALT have to be greater
#             if (refcount / sum) > 0.05 and (altcount / sum) > 0.05:
#                 curr_var1.append(record.CHROM)
#                 left1.append(record.POS)
#                 right1.append(record.POS + len(record.ALT))
#                 ref_seq1.append(record.REF)
#                 alt_seq1.append(str((record.ALT[0])))
#                 info.append(str(record.INFO))
#                 refct1.append(refcount)
#                 altct1.append(altcount)
#
#                 lists = curr_var1 + left1 + right1 + ref_seq1 + alt_seq1 + var_score + info
#                 normal_var.append(lists)
#
#
#             # throw away variant if not
#             else:
#                 pass
#
#
#     #transfer list into dataframe, it would be easier for following manipulation
#     normal=pd.DataFrame(normal_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score', 'info'])
#     tumor=pd.DataFrame(tumor_var,columns=['chrom','left','right','ref_seq','var_seq1','var_score', 'info'])
#
#     normal_seq2 = []
#     tumor_seq2 = []
#
#     normal.insert(5, "count1", altct1)
#     normal.insert(6, "count2", refct1)
#     for i in range(normal.shape[0]):
#
#         if normal['ref_seq'][i] > normal['var_seq1'][i]:
#             normal_seq2.append(normal['ref_seq'][i])
#
#         else:
#             normal_seq2.append(normal['var_seq1'][i])
#             normal['var_seq1'][i] = normal['ref_seq'][i]
#             # switching count values because the ref order is now var_seq1
#             temp = normal['count1'][i]
#             normal['count1'][i] = normal['count2'][i]
#             normal['count2'][i] = temp
#
#     normal.insert(5, "var_seq2", normal_seq2)  # inserts column with information given
#
#     tumor.insert(5, "count1", altct)
#     tumor.insert(6, "count2", refct)
#     for j in range(tumor.shape[0]):
#         if tumor['ref_seq'][j] > tumor['var_seq1'][j]:
#             tumor_seq2.append(tumor['ref_seq'][j])
#
#         else:
#             tumor_seq2.append(tumor['var_seq1'][j])
#             tumor['var_seq1'][j] = tumor['ref_seq'][j]
#             # switching count values because the ref order is now var_seq1
#             temp = tumor['count1'][j]
#             tumor['count1'][j] = tumor['count2'][j]
#             tumor['count2'][j] = temp
#
#     tumor.insert(5, "var_seq2", tumor_seq2)
#
#     normal = addALT_Seq(normal)
#     tumor = addALT_Seq(tumor)
#
#     #inserts and sets the var_index using a list of number within the range of the rows
#     tumor.insert(0, "var_index", list(range(tumor.shape[0])))
#     normal.insert(0, "var_index", list(range(normal.shape[0])))
#     #name = "Patient" + str(counts)
#     fil = inputfile.split(".")
#     id = fil[1] # Change index depending on the extracted file names (i.e. 0 for AML and ALL but 1 for THCA and thca)
#
#     ## inserts VCF ID as a column
#     tumor.insert(tumor.shape[1], "VCF_ID", id)
#     normal.insert(normal.shape[1], "VCF_ID", id)
#
#
#
#     ################################
#     # Inserts three columns, vcf file name, patient ID, and the alternate ID
#     ele1, ele2, ele3 = CaseIDs(inputfile)
#
#     tumor.insert(tumor.shape[1], "Patient_ID", ele3)
#     tumor.insert(tumor.shape[1], "Case_ID", ele2)
#
#     normal.insert(normal.shape[1], "Patient_ID", ele3)
#     normal.insert(normal.shape[1], "Case_ID", ele2)
#
#     ################################
#
#
#
#     ################################
#     # for VCF with info columns, use the previously defined functions to split this column
#     # if info column contains no information, comment that code out below
#     tumor = splitINFO(tumor)
#     tumor = parseCSQ(tumor, header)
#
#     normal = splitINFO(normal)
#     normal = parseCSQ(normal, header)
#
#     # Remove the Info and CSQ Column
#     tumor = tumor.drop(columns=["info", "CSQ"])
#     normal = normal.drop(columns=["info","CSQ"])
#     ##################################
#
#
#
#     ############### export the CSV file ###############
#     # tumor.to_csv(id + "_tumor.csv")
#     # normal.to_csv(id + "_normal.csv")
#     return tumor, normal
#
#
# ## perform vcf to csv, merged the many dataframes together without exporting the results into a csv
#
# import os
# import vcf
# import pandas as pd
# import re
# import itertools
# import numpy as np
# from numpy import log as ln
# import math
# from functools import reduce
# from statistics import mean
# import numpy as np
# import random
#
#
# # thca
# os.chdir("C:/Users/omodo/PycharmProjects/pythonProject/Fall 2024/Thyriod Cancer/Thyroid/Thyroid VCF files")
# files = os.listdir()
# thcaTumor = pd.DataFrame()
# thcaNormal = pd.DataFrame()
# thcaMerge = pd.DataFrame()
#
# for file in files:
#     thcaT, thcaN = writecsv(file, thcaheader)
#     thcaTumor = pd.concat([thcaTumor, thcaT], axis=0)
#     thcaNormal = pd.concat([thcaNormal, thcaN], axis=0)
#     thcaMerge = pd.concat([thcaMerge, thcaT], axis=0)
#     thcaMerge = pd.concat([thcaMerge, thcaN], axis=0)
#
#
# # Remove the chrM variants, those are to be ignored
# thcaTumor = thcaTumor[thcaTumor["chrom"] != "chrM"]
# thcaNormal = thcaNormal[thcaNormal["chrom"] != "chrM"]
#
#
#
# ## remove duplicates based on patient id and variant info
# thcaNormal = thcaNormal.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"])
# thcaTumor = thcaTumor.drop_duplicates(subset=["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"])
#
#
# ## subset the variant info and patients info only
# thcaptN = thcaNormal[["chrom", "left", "right", "ref_seq", "alt_seq", "VCF_ID", "Patient_ID"]]
# thcaptT = thcaTumor[["chrom", "left", "right", "ref_seq", "alt_seq", "VCF_ID", "Patient_ID"]]
# #
#
# ## group by variant info and insert N# and T# counts for each
# thcaptN = thcaptN.groupby(["Patient_ID"], as_index=False).agg(list)
#
# print("normal count")
# print(thcaptN)
#
# thcaptT = thcaptT.groupby(["Patient_ID"], as_index=False).agg(list)
#
# print("tumor count")
# print(thcaptT)
#
# thcaptN.insert(thcaptN.shape[1], "N#", 0)
# thcaptT.insert(thcaptT.shape[1], "T#", 0)
#
#
#
#
# def N_T(csv, col, outcol):
#
#     for row in range(csv.shape[0]):
#         csv[outcol][row] = len(csv[col][row])
#
#     return csv
#
#
# thcaptN = N_T(thcaptN, "chrom", "N#")
# thcaptT = N_T(thcaptT, "chrom", "T#")
#
#
# thcaptN = thcaNormal
# thcaptT = thcaTumor
#
# genomic_data = pd.merge(thcaptT, thcaptN, how="outer", on=["Patient_ID"])
# genomic_data.to_csv("Genomic Data.csv", index="False")
#
#
#
# ###Obtaining clinical data
# gdc_cancer = pd.read_csv("thca_tcga_clinical_data.csv")
# tcga_cancer = pd.read_csv("thpa_tcga_gdc_clinical_data.csv")
#
# clinical_files = ['thca_tcga_clinical_data.csv', 'thpa_tcga_gdc_clinical_data.csv']  # replace with actual file paths
# match_counts = []
#
# for clinical_file in clinical_files:
#     clinical_df = pd.read_csv(clinical_file)
#     clinical_df.columns = clinical_df.columns.str.strip()
#     merged_df = pd.merge(clinical_df, clinical_df, right_on='Patient ID',left_on="Patient_ID", how='inner')
#
#     match_counts[clinical_file] = len(merged_df)
#     print(f"Matches with {clinical_file}: {len(merged_df)}")
#
# # Find the clinical dataset with the most matches
# best_match = max(match_counts, key=match_counts.get)
# print(f"\nClinical dataset with most matches: {best_match} ({match_counts[best_match]} matches)")
#
# merged_df = pd.merge(gdc_cancer, tcga_cancer, right_on='Patient ID', left_on="Patient ID", how='inner')
# print(merged_df)
# merged_df = merged_df.loc[:, ~merged_df.columns.str.endswith('_y')]
#
# # Rename remaining columns by removing the '_x' suffix
# merged_df.columns = merged_df.columns.str.replace('_x', '', regex=False)
# clinical_data = merged_df[merged_df['Cancer Type Detailed'] == 'Papillary Thyroid Cancer']
# print(clinical_data)
# clinical_data = clinical_data.dropna(axis=1, how='all')
# print(clinical_data)
# clinical_data.to_csv("Clinical Dataset.csv", index=False)
#
# genomic_data = pd.read_csv("Thyroid_Patient_Check(in).csv")
# merge = pd.merge(clinical_data, genomic_data, how="inner", left_on=["Patient ID"],right_on=["Patient_ID"])

#
# #Check Mutation Count average
# combined_mean = merge[['N#', 'T#']].mean().mean()  # Mean across both columns
# combined_min = merge[['N#', 'T#']].stack().min()   # Minimum across both columns
# combined_max = merge[['N#', 'T#']].stack().max()   # Maximum across both columns
#
# print("Combined Mean:", combined_mean)
# print("Combined Minimum:", combined_min)
# print("Combined Maximum:", combined_max)
# stage_averages = merge.groupby('Neoplasm Disease Stage American Joint Committee on Cancer Code')[['N#', 'T#']].mean()
# print(stage_averages)

###Cox propotional hazard
# file_path = "Clinical data with mutation count.csv"
# data = pd.read_csv(file_path)
# data = data[[
#     "Overall Survival (Months)",  # Survival time
#     "Overall Survival Status",   # Event status (1 = death, 0 = censored)
#     "Diagnosis Age",             # Covariate: Age
#     "T#",            # Covariate: Total mutations
#     #"TMB (nonsynonymous)",       # Covariate: Tumor Mutational Burden
#     #"Prior Cancer Diagnosis Occurrence",
#     #"Prior Cancer Diagnosis Occurrence ",
#     #"Prior Malignancy",
#     #"Prior Treatment",
#     "American Joint Committee on Cancer Tumor Stage Code",  # Covariate: Tumor stage
#     "Neoplasm Disease Stage American Joint Committee on Cancer Code",  # Covariate: Disease stage
#     "Sex",                       # Covariate: Gender
#     "Race Category"  ,
#     "Ethnicity Category"# Covariate: Ethnicity,
#     "DiseaseFree Status "
# ]]

# output_path = "subset_data.csv"
# subset_data.to_csv(output_path, index=False)
# print(f"Subset data saved to {output_path}")

# # Define mapping dictionaries
# stage_mapping = {"Stage IVA": 5, "Stage II": 2, "Stage I": 1, "Stage III": 3, "Stage IV": 4}
# version_mapping = {
#     "T1": 1, "T1a": 1.1, "T1b": 1.2,  # Subcategories of T1
#     "T2": 2,
#     "T3": 3,
#     "T4": 4, "T4a": 4.1, "T4b": 4.2,  # Subcategories of T4
#     "TX": 0  # Assuming TX is unknown, map to 0
# }
#
# print(data["American Joint Committee on Cancer Tumor Stage Code"].unique())
# print(data["Neoplasm Disease Stage American Joint Committee on Cancer Code"].unique())
#
# # Apply mapping to columns
# data["American Joint Committee on Cancer Tumor Stage Code"] = data["American Joint Committee on Cancer Tumor Stage Code"].map(version_mapping)
# data["Neoplasm Disease Stage American Joint Committee on Cancer Code"] = data["Neoplasm Disease Stage American Joint Committee on Cancer Code"].map(stage_mapping).fillna(-1)
#
# #"Prior Malignancy","Prior Cancer Diagnosis Occurrence ","Prior Treatment",
# # Convert categorical variables to dummy/indicator variables
# data = pd.get_dummies(data, columns=["Ethnicity Category", "Sex", "Race Category"], drop_first=True)
# data = data.drop(columns=["Race Category_ASIAN"])
#
# print(data)
# cph = CoxPHFitter()
#
# # Fit the Cox model to the data
# cph.fit(data, duration_col='Overall Survival (Months)', event_col='Overall Survival Status')
# cph.print_summary()
# cph.summary.to_csv('cox_model_results.csv')


# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from lifelines import CoxPHFitter
#
# # Load datasets
# snv_data = pd.read_csv("snv_data.csv")  # Columns: Sample_ID, Gene, Impact, etc.
# clinical_data = pd.read_csv("clinical_data.csv")  # Columns: Sample_ID, Survival_Time, Event
#
# # Parameters
# frequency_threshold = 5  # Define mutation frequency threshold
# impact_filter = ['HIGH', 'MODERATE']  # Functional impact categories to include
#
# # Step 1: Mutation Frequency Analysis
# mutation_counts = snv_data['Gene'].value_counts()
# high_freq_mutations = mutation_counts[mutation_counts > frequency_threshold]
# high_freq_genes = high_freq_mutations.index.tolist()
#
# # Step 2: Functional Effect Filtering
# filtered_mutations = snv_data[snv_data['Impact'].isin(impact_filter)]
# combined_candidates = filtered_mutations[filtered_mutations['Gene'].isin(high_freq_genes)]
#
# # Step 3: Prioritize Mutations
# # Assign scores: Frequency x Impact (numerical value for impact categories)
# impact_scores = {'HIGH': 2, 'MODERATE': 1}  # Scoring for impacts
# combined_candidates['Frequency'] = combined_candidates['Gene'].map(mutation_counts)
# combined_candidates['Impact_Score'] = combined_candidates['Impact'].map(impact_scores)
# combined_candidates['Score'] = (
#     combined_candidates['Frequency'] * combined_candidates['Impact_Score']
# )
# prioritized_mutations = combined_candidates.sort_values(by='Score', ascending=False)
#
# # Step 4: Visualization
# plt.figure(figsize=(12, 6))
# sns.scatterplot(
#     data=prioritized_mutations,
#     x='Frequency',
#     y='Impact_Score',
#     size='Score',
#     hue='Gene',
#     alpha=0.7
# )
# plt.title('Mutation Frequency vs Functional Impact')
# plt.xlabel('Frequency')
# plt.ylabel('Impact Score')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.tight_layout()
# plt.show()
#
# # Step 5: Link Mutations to Clinical Data
# for gene in prioritized_mutations['Gene'].unique():
#     clinical_data[f'{gene}_mutation'] = clinical_data['Sample_ID'].apply(
#         lambda x: 1 if x in combined_candidates[combined_candidates['Gene'] == gene]['Sample_ID'].unique() else 0
#     )
#
# # Step 6: Survival Analysis
# cph = CoxPHFitter()
# cph.fit(clinical_data, duration_col='Survival_Time', event_col='Event')
# cph.print_summary()  # Print Cox model results
#
# # Save results
# prioritized_mutations.to_csv("prioritized_mutations.csv", index=False)
# clinical_data.to_csv("updated_clinical_data.csv", index=False)
#
#
# file_path = "Clinical data with mutation count.csv"  # Replace with the actual file path
# data = pd.read_csv(file_path)
# #data = data[["Patient ID","American Joint Committee on Cancer Tumor Stage Code","Sex","Race Category","Diagnosis Age"]]
#
# thyriod = pd.read_csv("C:/Users/omodo/PycharmProjects/pythonProject/Fall 2024/Thyriod Cancer/thca_check_cut_off.csv")
# thyriod.drop_duplicates().reset_index()
# thyriod = thyriod[["chrom","left","ref_seq","alt_seq","right_x","Info_SIFT_Score","Info_PolyPhen_Score","FATHMM_Score","Amino_acids","Patient_ID_Normal","Patient_ID_Tumor","Patient_ID","SYMBOL"]]
# #thyriod = thyriod[["chrom","left","ref_seq","alt_seq","Patient_ID","SYMBOL"]]
#
# check = pd.merge(data,thyriod, how="left", left_on="Patient ID", right_on="Patient_ID")
# check = check.groupby(["SYMBOL"], as_index=False).agg(list)
#
# check['Patient Count'] = check["Patient ID"].apply(lambda x: len(set(x)))
# check['Mutation Count'] = check["chrom"].apply(len)
# check.to_csv("Check mutaion pt count.csv", index=False)
#
#
# summary = data.groupby('Neoplasm Disease Stage American Joint Committee on Cancer Code').agg(
#     avg_normal_mutation=('N#', 'mean'),
#     sd_normal_mutation=('N#', 'std'),
#     avg_tumor_mutation=('T#', 'mean'),
#     sd_tumor_mutation=('T#', 'std')
# ).reset_index()
#
#
# summary['normal_mutation_summary'] = (
#     summary['avg_normal_mutation'].round(2).astype(str) + " ± " + summary['sd_normal_mutation'].round(2).astype(str)
# )
# summary['tumor_mutation_summary'] = (
#     summary['avg_tumor_mutation'].round(2).astype(str) + " ± " + summary['sd_tumor_mutation'].round(2).astype(str)
# )
#
# # Display the final summary
# print(summary[['Neoplasm Disease Stage American Joint Committee on Cancer Code', 'normal_mutation_summary', 'tumor_mutation_summary']])
#
# #Display the results
# print(summary)
#
# # Optionally save to a new CSV file
# summary.to_csv("mutation_summary.csv", index=False)

# For TMB burden

file_path = "Clinical data with mutation count.csv"  # Replace with the actual file path
data = pd.read_csv(file_path)

summary = data.groupby('Neoplasm Disease Stage American Joint Committee on Cancer Code').agg(
    avg_normal_mutation=('TMB (nonsynonymous)', 'mean'),
    sd_normal_mutation=('TMB (nonsynonymous)', 'std')
).reset_index()


summary['normal_mutation_summary'] = (
    summary['avg_normal_mutation'].round(2).astype(str) + " ± " + summary['sd_normal_mutation'].round(2).astype(str)
)

# Display the final summary
print(summary[['Neoplasm Disease Stage American Joint Committee on Cancer Code', 'normal_mutation_summary']])

#Display the results
print(summary)

# Optionally save to a new CSV file
summary.to_csv("mutation_summary.csv", index=False)





# #thyriod = thyriod[["Patient_ID","SYMBOL",'FATHMM_Score',"Info_PolyPhen_Score","Info_SIFT_Score"]]
# check = pd.merge(data,thyriod, how="left", left_on="Patient ID", right_on="Patient_ID")
#
#
# def get_top_100_genes(df, score_column, symbol_column):
#
#     sorted_df = df.sort_values(by=score_column, ascending=False)
#
#     unique_top = sorted_df.drop_duplicates(subset=symbol_column)
#
#     return unique_top
#
#
# # Specify the relevant columns
# fathmm_score_column = 'FATHMM_Score'
# polyphen_score_column = 'Info_PolyPhen_Score'
# gene_symbol_column = 'SYMBOL'  # Replace with the actual column name for gene symbols
#
# # Get top 100 unique genes for FATHMM and PolyPhen
# fathmm_top_100 = get_top_100_genes(check, fathmm_score_column, gene_symbol_column)
# polyphen_top_100 = get_top_100_genes(check, polyphen_score_column, gene_symbol_column)
#
# # Save the results to separate files or continue processing
# fathmm_top_100.to_csv('fathmm_top_100_genes.csv', index=False)
#polyphen_top_100.to_csv('polyphen_top_100_genes.csv', index=False)
#

### code for patient count
#thyriod = thyriod[["chrom","left","ref_seq","alt_seq","right_x","Info_SIFT_Score","Info_PolyPhen_Score","FATHMM_Score","Amino_acids","Patient_ID_Normal","Patient_ID_Tumor","Patient_ID","SYMBOL"]]
# thyriod = thyriod[["chrom","left","ref_seq","alt_seq","Patient_ID","SYMBOL"]]
#
# check = pd.merge(data,thyriod, how="left", left_on="Patient ID", right_on="Patient_ID")


#
# thyriod1 = thyriod.groupby(["SYMBOL"], as_index=False).agg(list)
# thyriod1['Gene_Count'] = thyriod1["ref_seq"].apply(len)
# thyriod1['Patient Count'] = thyriod1["Patient_ID_Tumor"].apply(lambda x: len(set(x)))
#
#
# # Normalize mutation count by patient count
# thyriod1['Normalized_Mutation_Rate'] = thyriod1['Gene_Count'] / thyriod1['Patient Count']
#



#thyriod1.to_csv("gene.csv",index=False)

###Survival analysis
#
# # Filter patients with MUC4 mutation
# check['MUC4 Mutation'] = check['SYMBOL'].apply(lambda x: "Yes" if x == "MUC4" else "No")
#
# check = check.drop_duplicates()
# check = check.reset_index()
#
# check.to_csv("Check mutation.csv", index=False)
#
#
# # Calculate survival rates for MUC4 vs. non-MUC4 patients
# survival_by_muc4 = check.groupby('MUC4 Mutation')['Overall Survival Status'].value_counts(normalize=True).unstack()
#
# print(survival_by_muc4)
#
#


#KaplanMeierFitter for gene and patient survival
# from lifelines import KaplanMeierFitter
# import matplotlib.pyplot as plt
#
# # Split data by MUC4 mutation presence
# muc4_yes = check[check['MUC4 Mutation'] == 'Yes']
# muc4_no = check[check['MUC4 Mutation'] == 'No']
#
# # Kaplan-Meier analysis
# kmf = KaplanMeierFitter()
#
# # MUC4 mutation present
# kmf.fit(muc4_yes['Overall Survival (Months)'], event_observed=muc4_yes['Overall Survival Status'] == '1')
# ax = kmf.plot(label='MUC4 Mutation')
#
# # MUC4 mutation absent
# kmf.fit(muc4_no['Overall Survival (Months)'], event_observed=muc4_no['Overall Survival Status'] == '1')
# kmf.plot(ax=ax, label='No MUC4 Mutation')
#
# plt.title('Kaplan-Meier Survival Curve by MUC4 Mutation')
# plt.xlabel('Time (days)')
# plt.ylabel('Survival Probability')
# plt.legend()
# plt.show()
#
# # Group by tumor stage and MUC4 mutation status
# muc4_stage_distribution = check.groupby(['Neoplasm Disease Stage American Joint Committee on Cancer Code', 'MUC4 Mutation']).size().unstack(fill_value=0)
# print(muc4_stage_distribution)
# from lifelines import KaplanMeierFitter
# import matplotlib.pyplot as plt
#
# # Initialize the Kaplan-Meier estimator
# kmf = KaplanMeierFitter()
#
# # Plot survival curves for each stage and MUC4 mutation status
# plt.figure(figsize=(14, 10))
#
# # Iterate over each stage
# for stage in check['Neoplasm Disease Stage American Joint Committee on Cancer Code'].unique():
#     stage_data = check[check['Neoplasm Disease Stage American Joint Committee on Cancer Code'] == stage]
#
#     # Plot for MUC4 mutation present
#     muc4_patients_stage = stage_data[stage_data['MUC4 Mutation'] == 'Yes']
#     kmf.fit(muc4_patients_stage['Overall Survival (Months)'],
#             event_observed=muc4_patients_stage['Overall Survival Status'] == 1)
#     kmf.plot(label=f'Stage {stage} - MUC4 Mutation')
#
#     # Plot for MUC4 mutation absent
#     non_muc4_patients_stage = stage_data[stage_data['MUC4 Mutation'] == 'No']
#     kmf.fit(non_muc4_patients_stage['Overall Survival (Months)'],
#             event_observed=non_muc4_patients_stage['Overall Survival Status'] == 1)
#     kmf.plot(label=f'Stage {stage} - No MUC4 Mutation')
#
# plt.title('Survival Analysis by Tumor Stage and MUC4 Mutation Status')
# plt.xlabel('Time (Months)')
# plt.ylabel('Survival Probability')
# plt.legend()
# plt.show()#Tumor mutation burden
