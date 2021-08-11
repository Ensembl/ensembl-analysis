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
    df_genmark_final = pd.DataFrame(df.columns)
    for i in gene_list:

        df_subset = df_genmark[df_genmark["gene_id"] == str(i)]
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

        if i == list(gene_list)[0]:
            df_genmark_final = pd.concat(
                [df_subset.iloc[index - 1 :], line, df_subset.iloc[: index - 1]]
            ).reset_index(drop=True)
        else:
            df_genmark_final = df_genmark_final.append(
                pd.concat(
                    [df_subset.iloc[index - 1 :], line, df_subset.iloc[: index - 1]]
                ).reset_index(drop=True)
            ).reset_index(drop=True)

    df_genmark_final = df_genmark_final.replace(np.nan, "")
    
    return df_genmark_final


def analyse_augustus(df):
    df_augustus = df[df["source"] == "AUGUSTUS"].reset_index()
    df_augustus = df_augustus.drop(["index"], axis=1)
    df_augustus["score"] = df_augustus["score"].replace(np.nan, ".")
    df_augustus["frame"][df_augustus["feature"] == "exon"] = "."
    df_augustus["frame"][df_augustus["feature"] == "intron"] = "."
    df_augustus["frame"][df_augustus["feature"] == "gene"] = "."
    df_augustus["frame"][df_augustus["feature"] == "transcript"] = "."
    
    df_augustus = df_augustus.drop(
        df_augustus[df_augustus.duplicated()].index, axis=0
    ).reset_index(drop=True)
    df_augustus = df_augustus.drop(
        df_augustus[df_augustus["feature"] == "gene"].index, axis=0
    ).reset_index(drop=True)
    df_augustus = df_augustus.drop(
        df_augustus[df_augustus["feature"] == "intron"].index, axis=0
    ).reset_index(drop=True)
    df_augustus = df_augustus.drop_duplicates().reset_index(drop=True)
    df_augustus["transcript_id"] = [
        re.sub(r"^jg", "braker_augustus", str(x)) for x in df_augustus["transcript_id"]
    ]
    df_augustus["gene_id"] = [
        re.sub(r"^jg", "braker_augustus", str(x)) for x in df_augustus["gene_id"]
    ]
    transcript_list = set(df_augustus["transcript_id"])
    df_augustus_final = pd.DataFrame()
    for i in transcript_list:
        
        df_subset = df_augustus[df_augustus["transcript_id"] == i]

        old_index = df_subset[df_subset["feature"] == "transcript"].index.item()

        df_subset_new = df_subset.T.reindex(
            columns=(
                [old_index] + list([a for a in df_subset.T.columns if a != old_index])
            )
        )

        df_augustus_final = df_augustus_final.append(df_subset_new.T).reset_index(
            drop=True
        )
        

    
    return df_augustus_final


def create_file(df):
    new_file = []

    for l in range(0, len(df)):
        temp = []
        first_line_part = "\t".join([str(i) for i in df.loc[l, "seqname":"frame"]])
        second_line_part = df.loc[l, "transcript_id":].tolist()
        col_index = 8
        for t in second_line_part:
            if len(t) != 0:
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

    output = genmark_list+augustus_list
    save_file(output, output_file)
