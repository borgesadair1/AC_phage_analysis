#!/usr/bin/env python
import subprocess
import pandas as pd
from Bio import SeqIO
import argparse


def eval_AC(cd4, cd11, cd15,Len_kb):

    answer_list = []
    counter = 0
    
    
    for val in cd11:
        answer = 11 
        code4_diff = cd4[counter] - val 
        code15_diff = cd15[counter] - val
        len_kb = Len_kb[counter]
        
        if len_kb >= 100: 
            cutoff = 0.05
        if len_kb < 100:
            cutoff = 0.1
            
        if (code4_diff > code15_diff) & (code4_diff > cutoff):
            answer = 4
        
        if (code15_diff > code4_diff) & (code15_diff > cutoff):
            answer = 15

        counter +=1
 
        answer_list.append(int(answer))
  

    
    return(answer_list)

def main():
    parser = argparse.ArgumentParser(description='''calculates coding density after prodigal gene prediction ''')
    parser.add_argument('-sc','--scaffold', help='Path to scaffold fasta file.',required=True) 
    parser.add_argument('-o', '--output_directory', help='Path to output directory.', required=True)
    parser.add_argument('-c4', '--code4', help='Path to faa file predicted with gcode4', required=True)
    parser.add_argument('-c15', '--code15', help='Path to faa file predicted with gcode15', required=True)
    parser.add_argument('-c11', '--code11', help='Path to faa file predicted with gcode11', required=True)

    
    args = parser.parse_args()
    in_fasta = args.scaffold
    out_dir = args.output_directory
    out_faa_4 = args.code4
    out_faa_11 = args.code11
    out_faa_15 = args.code15
    
    file_dict = {out_faa_11: 'code11', out_faa_15: 'code15', out_faa_4: 'code4'}
    aa_dict = {}
    
    for file in file_dict.keys():
        
        for record in SeqIO.parse(file, "fasta"):
            code = file_dict[file]
            aa = len(record.seq)
            gene_name = record.description.split(" ")[0]
            scaff_name = gene_name.split("_")[:-1]
            scaff_name = "_".join(scaff_name) + "." + code
            
            if scaff_name in aa_dict.keys():
                aa_dict[scaff_name] = aa_dict[scaff_name] + aa
            else:
                aa_dict[scaff_name] = aa
    aa_df = pd.DataFrame.from_dict(aa_dict, orient='index',).reset_index()
    aa_df = aa_df.rename({'index': 'Scaffold', 0: 'Len_AA'}, axis=1)

    aa_df["Code"] = aa_df["Scaffold"].apply(lambda x: x.split(".code")[1])
    aa_df["Scaffold"] = aa_df["Scaffold"].apply(lambda x: x.split(".code")[0])
    aa_df["Len_AA"] = aa_df["Len_AA"].apply(lambda x: x*3)


    scaff_dict = {}
        
    for record in SeqIO.parse(in_fasta, "fasta"):
        scaff_name = record.description.split(" ")[0]
        nt = len(record.seq)
        scaff_dict[scaff_name] = nt



    scaff_df = pd.DataFrame.from_dict(scaff_dict, orient='index',).reset_index()
    
    scaff_df = scaff_df.rename({'index': 'Scaffold', 0: 'Len_nt'}, axis=1)
    merged = aa_df.merge(scaff_df, left_on = "Scaffold", right_on = "Scaffold", how = "left")
    merged['CD'] = merged['Len_AA']/merged['Len_nt']
	
    code4_df = merged.loc[merged["Code"] == "4"]
    code11_df = merged.loc[merged["Code"] == "11"]
    code15_df = merged.loc[merged["Code"] == "15"]
	
    code4_df = code4_df.rename({'CD': 'Code4_CD'}, axis=1) 
    code11_df = code11_df.rename({'CD': 'Code11_CD'}, axis=1) 
    code15_df = code15_df.rename({'CD': 'Code15_CD'}, axis=1) 
	
    merged_codes_df = code4_df.merge(code11_df, on = "Scaffold")
    merged_codes_df = merged_codes_df.merge(code15_df, on = "Scaffold")
    merged_codes_df = merged_codes_df[["Scaffold", "Len_nt_x", "Code4_CD", "Code11_CD", "Code15_CD"]]
	
    merged_codes_df = merged_codes_df.rename({'Len_nt_x': 'Len_nt'}, axis=1)
    merged_codes_df['Len_kb'] = merged_codes_df['Len_nt']/1000
	
    merged_codes_df["pred_code"] = eval_AC(merged_codes_df["Code4_CD"].to_list(),merged_codes_df["Code11_CD"].to_list(),merged_codes_df["Code15_CD"].to_list(), merged_codes_df["Len_kb"].to_list()) 

    merged_codes_df.to_csv(out_dir + "/CD_analysis.csv", index = False)







if __name__ == "__main__":
    main()
