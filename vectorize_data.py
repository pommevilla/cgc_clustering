#!/usr/bin/env python
# ---------------------------
# Vectorizes the PULs in the xls by their enzymes
# This creates two versions of the vectorized data:
#  1. One where the domains are split by comma
#  2. One where the domains are kept together
# This is because there might be some biological significance
# to the domains being kept together
# Author: Paul Villanueva (github.com/pommevilla)
# ---------------------------
import pandas as pd


def vectorize_df(this_df):
    """
    Vectorizes a PUL dataframe by the enzymes.
    """
    exploded = this_df.explode("enzyme_list")
    exploded = exploded.rename(columns={"enzyme_list": "enzyme"})

    exploded["presence"] = 1

    vectorized_cgcs = exploded.pivot_table(
        index="pul_id",
        columns="enzyme",
        values="presence",
        margins=True,
        aggfunc="sum",
        margins_name="total_enzymes",
        fill_value=0,
    ).drop("total_enzymes")

    vectorized_cgcs.columns = vectorized_cgcs.columns.get_level_values(0)

    big_df = vectorized_cgcs.join(this_df.set_index("pul_id"))

    return big_df


if __name__ == "__main__":
    df = pd.read_excel(
        "data/dbCAN-PUL.substrate.mapping.xls",
        sheet_name="Table S1",
        skiprows=1,
        names=[
            "pul_id",
            "pmid",
            "substrate_04142023",
            "substrate_09012022",
            "substrate_07012022",
            "predicted_enzymes",
            "notes",
        ],
    )

    # This df splits comma-separated domains
    less_domains = df.copy()
    less_domains["enzyme_list"] = (
        less_domains["predicted_enzymes"].str.replace(",", "|").str.split("|")
    )
    less_domains["enzyme_count"] = less_domains["enzyme_list"].str.len()

    # This df keeps comma-separated domains together
    more_domains = df.copy()
    more_domains["enzyme_list"] = more_domains["predicted_enzymes"].str.split("|")
    more_domains["enzyme_count"] = more_domains["enzyme_list"].str.len()

    filenames = ["less_domains.csv", "more_domains.csv"]
    dfs = [less_domains, more_domains]

    for fout, dataframe in zip(filenames, dfs):
        vectorized_df = vectorize_df(dataframe)
        print(f"Dims of {fout} after vectoring: {dataframe.shape}")
        print(dataframe.head())
        vectorized_df.to_csv(f"output/vectorized_puls/{fout}")
