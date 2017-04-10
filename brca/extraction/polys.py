import pandas as pd
import re
from get_attributes import get_genome_hgvs, format_p_hgvs
import time
import os


def get_full_protein_hgvs(shortened_hgvs):
    cwd = os.getcwd()
    """Takes input of e.g. M123A and returns full nomenclature e.g. p.Met123Ala.

     Reads text file of amino acid codes into dataframe and separates into two lists: single-letter codes and
    three-letter codes and uses these to create a reference dictionary.
     Converts the input into a list for easier modification (e.g. ['M', '1', '2', '3', 'A']) and checks each item for a
    matching key in the reference dictionary. 'M' and 'A' would match to 'Met' and 'Ala', so these replace the single
    letters. Then a "p." is added to the start of the list, and the list is converted back to a string.
     If the second match is the same as the first (i.e. amino acids are the same) the second is converted to "="
    instead, to comply with HGVS formatting guidelines.

    :param shortened_hgvs:e.g. M123A
    :return: full_hgvs: e.g. p.Met123Ala
    """
    aa_df = pd.read_table("{cwd}/brca/extraction/amino_acids.txt".format(cwd=cwd), header=None, names=['Single', 'Triple'], sep='=')
    single = aa_df['Single'].tolist()
    triple = aa_df['Triple'].tolist()
    aa_dict = dict(zip(single, triple))
    syn_list = []

    hgvs_list = list(shortened_hgvs)

    for item in hgvs_list:
        item_index = hgvs_list.index(item)
        if item.upper() in aa_dict:
            if item in syn_list:
                aa = "="
            else:
                syn_list.append(item)
                aa = aa_dict[item.upper()]
            hgvs_list[item_index] = aa
    hgvs_list[:0] = ["p."]
    full_hgvs = format_p_hgvs("".join(hgvs_list))

    return full_hgvs


def get_brca_polys(variant_file):
    """Gets gene, cdna and protein nomenclature for all class 1 and 2 variants listed in Excel spreadsheet.

    Works with layout of WMRGL BRCA mutation list:
        A: Gene - row[0] in poly_df dataframe
        D: HGVS Nucleotide - row[1] in poly_df dataframe
        E: Codon - row[2] in poly_df dataframe

    (1) Captures sheet name with "Poly" in title for use in read_excel (exception raised if none found).
    (2) Read the three relevant columns into a dataframe using pandas.
    (3) For each row in the dataframe, extract the gene and cDNA and protein changes and write into a CSV file.
            (i) Check each row is valid, if not raise an exception.
            (ii) Some cDNA entries are written in two ways so split these and extract only the first one.
            (iii) If a protein change is provided, convert to full HGVS nomenclature.
            (iv) Extract the gene.

    :param variant_file: Excel document with list of variants, sheet with polys should have "Poly" in sheet name.
    """
    sheetname = None
    special_characters = [";", ":", "!", "?", ","]
    xl = pd.ExcelFile(variant_file)
    sheet_names = xl.sheet_names
    for name in sheet_names:
        if "Poly" in name:
            sheetname = name
        else:
            pass
    if sheetname:
        poly_df = pd.read_excel(
            variant_file,
            header=0,
            parse_cols="A, D, E",
            skiprows=2,
            names=["Gene", "cDNA", "Protein"],
            sheetname=sheetname,
            index_col=None
        )
        poly_df = poly_df.where((pd.notnull(poly_df)), None)  # change NaN to None, easier to work with.
        poly_df = poly_df.dropna(thresh=2)
        final_df = pd.DataFrame(
            columns=["Gene", "RefSeq Transcript ID", "Genomic Change", "cDNA Change", "Protein Impact"]
        )
        time.sleep(0.01)
        for row_index, row in poly_df.iterrows():
            for item in special_characters:
                if row[0] and item in row[0]:
                    raise Exception(
                        "Unexpected {character} in row {index}. Please remove and try again.".format(
                            character=item, index=row_index + 4
                        )
                    )
                if row[1] and item in row[1]:
                    raise Exception(
                        "Unexpected {character} in row {index}. Please remove and try again.".format(
                            character=item, index=row_index + 4
                        )
                    )
                if row[2] and item in row[2]:
                    raise Exception(
                        "Unexpected {character} in row {index}. Please remove and try again.".format(
                            character=item, index=row_index + 4
                        )
                    )
            if row[2] and row[1] is None:
                raise Exception(
                    "Entry {row_index} has a protein HGVS recorded but no cDNA HGVS. "
                    "Please fix and try again.".format(row_index=row_index+4)
                )
            else:
                split_terms = row[1].split()
                cdna = split_terms[0].strip("()")
            if row[1] and row[2] is None:
                protein = ""
            elif "p." in row[2]:
                protein = row[2]
            elif re.match(r"[a-zA-Z]{3}", row[2]):
                protein = "p." + row[2]
            else:
                split_terms = row[2].split(' ')
                protein = get_full_protein_hgvs(split_terms[0])
            gene = row[0].replace(" ", "")
            if gene == "BRCA1":
                refseq = "NM_007294.3"
            elif gene == "BRCA2":
                refseq = "NM_000059.3"
            g_hgvs = get_genome_hgvs("hg19", refseq, cdna)
            var_df = pd.DataFrame(
                [[gene, refseq, g_hgvs, cdna, protein]],
                columns=["Gene", "RefSeq Transcript ID", "Genomic Change", "cDNA Change", "Protein Impact"]
            )
            final_df = final_df.append(var_df)
    else:
        raise Exception("Keyword \"Poly\" not found in any sheet names in Excel document. Please fix and try again.")
    return final_df


