from gtfparse import read_gtf
import numpy as np
import pandas as pd
import sys
import logging
import re

logging.basicConfig(
    filename="app.log",
    filemode="w",
    format="%(asctime)s,%(msecs)d %(name)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
    level=logging.DEBUG,
)


def analyse_genmark(df):
    df_genmark = df[df["source"] == "GeneMark.hmm"].reset_index()
    df_genmark = df_genmark.drop(["index"], axis=1)
    df_genmark["score"] = df_genmark["score"].replace(np.nan, ".")
    df_genmark["frame"][df_genmark["feature"] == "exon"] = "."
    # dropping ALL duplicate values
    df_genmark = df_genmark.drop_duplicates().reset_index(drop=True)
    df_genmark["transcript_id"] = [
        re.sub(r"[_t]\b", "braker_transcript_GenMark", str(x))
        for x in df_genmark["transcript_id"]
    ]
    df_genmark["gene_id"] = [
        re.sub(r"[_g]\b", "braker_gene_GenMark", str(x)) for x in df_genmark["gene_id"]
    ]
    
    gene_list = set(df_genmark["gene_id"])
    df_genmark_final = pd.DataFrame()
    for i in gene_list:

        df_subset = df_genmark[df_genmark["gene_id"] == str(i)]
        if len(df_subset)>1:
            start = min(df_subset["start"])
            stop = max(df_subset["end"])
            index = min(df_subset.index)
            line = pd.DataFrame(
                {
                "seqname": df_subset["seqname"][index],
                "source": df_subset["source"][index],
                "feature": "transcript",
                "start": start,
                "end": stop,
                "score": ".",
                "strand": df_subset["strand"][index],
                "frame": ".",
                "transcript_id": df_subset["transcript_id"][index],
                "gene_id": df_subset["gene_id"][index],
                },
                index=[index],
            )
            df_subset = pd.concat([df_subset.iloc[index - 1 :], line, df_subset.iloc[: index - 1]]).reset_index(drop=True)                                                  
            df_subset.fillna('',inplace=True)    
            df_subset["score"] = df_subset["score"].replace(0.0, ".")
            old_index = df_subset[df_subset["feature"] == "transcript"].index.item() 
            df_subset_new = df_subset.T.reindex(columns=([old_index] + list([a for a in df_subset.T.columns if a != old_index])
                        )
                )
            df_genmark_final = df_genmark_final.append(df_subset_new.T).reset_index(
                    drop=True
                )
    
    return df_genmark_final


def analyse_augustus(df):
    df_augustus = df[df["source"] == "AUGUSTUS"].reset_index()
    df_augustus = df_augustus.drop(["index"], axis=1)
    df_augustus["score"] = df_augustus["score"].replace(np.nan, ".")
    df_augustus["frame"][df_augustus["feature"] == "exon"] = "."
    df_augustus["frame"][df_augustus["feature"] == "intron"] = "."
    df_augustus["frame"][df_augustus["feature"] == "gene"] = "."
    df_augustus["frame"][df_augustus["feature"] == "transcript"] = "."
    df_augustus["transcript_id"] = [
        re.sub(r"^g", "braker_augustus", str(x)) for x in df_augustus["transcript_id"]
    ]
    df_augustus["gene_id"] = [
        re.sub(r"^g", "braker_augustus", str(x)) for x in df_augustus["gene_id"]
    ]    

    df_augustus = df_augustus.drop(df_augustus[df_augustus["feature"] == "gene"].index, axis=0).reset_index(drop=True)
    df_augustus = df_augustus.drop(df_augustus[df_augustus["feature"] == "intron"].index, axis=0).reset_index(drop=True)

    for i in df_augustus[df_augustus['feature']=='transcript'].index.tolist():
        df_augustus['transcript_id'][i]=df_augustus['transcript_id'][i-1]
        df_augustus['gene_id'][i]=df_augustus['gene_id'][i-1]

    df_augustus_final = pd.DataFrame()

    for i in set(df_augustus["gene_id"]):
        df_gene = df_augustus[df_augustus["gene_id"] == i]
        for j in set(df_gene["transcript_id"]):

                
            df_transcript = df_gene[df_gene["transcript_id"] == j]

            if len(df_transcript[df_transcript["feature"] == "transcript"]) > 1 or len(df_transcript[df_transcript["feature"] == "transcript"]) ==0:

                df_transcript=df_transcript.drop(df_transcript[df_transcript["feature"] == "transcript"].index, axis=0).reset_index(drop=True)           
                start=min(df_transcript['start'])
                stop=max(df_transcript['end'])
                index=min(df_transcript.index)          
                line = pd.DataFrame({"seqname":df_transcript['seqname'][index],"source":df_transcript['source'][index], "feature":'transcript',"start":start, "end":stop,"score":'.', "strand":df_transcript['strand'][index], "frame":'.', "transcript_id":df_transcript['transcript_id'][index], "gene_id":df_transcript['gene_id'][index],}, index=[index])                
                df_transcript=pd.concat([df_transcript.iloc[index-1:], line, df_transcript.iloc[:index-1]]).reset_index(drop=True)
                df_transcript.fillna('',inplace=True)
            print(df_transcript)
            old_index = df_transcript[df_transcript["feature"] == "transcript"].index.item()
            df_subset_new = df_transcript.T.reindex(
                columns=(
                    [old_index] + list([a for a in df_transcript.T.columns if a != old_index])
                        )
                )

            df_augustus_final = df_augustus_final.append(df_subset_new.T).reset_index(
                    drop=True
                )

    return df_augustus_final

def create_file(df):
    new_file = []
    df.fillna('',inplace=True)
    for l in range(0, len(df)):
        temp = []
        first_line_part = "\t".join([str(i) for i in df.loc[l, "seqname":"frame"]])
        second_line_part = df.loc[l, "transcript_id":].tolist()
        col_index = 8
        #print(second_line_part)
        for t in second_line_part:
            #print(t)
            if t!='':
            #if len(t) != 0:
                temp.append(df.columns.values.tolist()[col_index] + " " + '"' + t + '"')
                col_index = col_index + 1

        temp_columns = "; ".join([str(i) for i in temp])
        final = first_line_part + "\t" + temp_columns + ";\n"

        new_file.append(final)
    return new_file


def save_file(output, output_file):
    tfile = open(output_file, "a")
    for i in output:
        tfile.write(i)

    tfile.close()


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    logging.info(input_file)
    df = read_gtf(str(input_file))
    source = set(df["source"])

    if "GeneMark.hmm" in source:
        
        df_genmark_final = analyse_genmark(df)
        genmark_list = create_file(df_genmark_final)

    if "AUGUSTUS" in source:
        df_augustus_final = analyse_augustus(df)
        augustus_list = create_file(df_augustus_final)

    output = genmark_list + augustus_list
    save_file(output, output_file)
