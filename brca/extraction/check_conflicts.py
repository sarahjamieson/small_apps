import pandas as pd


def check_conflicts(dataframe, error_file):
    """Checks a dataframe and identifies variants which have been classified and reported differently.

    :param error_file: txt file for writing the conflicts to.
    :param dataframe: should have columns labelled "cDNA Change", "Gene" and "Assigned Pathogencity Score".
    """
    var_list = []
    class_list = []
    repeats = 0
    conflicts = 0

    with open(error_file, "ab+") as errors:
        for row_index, row in dataframe.iterrows():
            if row["cDNA Change"] and not pd.isnull(row['cDNA Change']):
                lookup = "{gene}_{cdna}".format(gene=row["Gene"], cdna=row['cDNA Change'])
                if lookup in var_list:
                    repeats += 1
                    var_index = var_list.index(lookup)
                    classification = class_list[var_index]
                    if classification != row['Assigned Pathogenicity Score']:
                        conflicts += 1
                        errors.write("{cdna} in {gene} reported as {score1} and {score2}.\n".format(
                            cdna=row["cDNA Change"], gene=row["Gene"], score1=row['Assigned Pathogenicity Score'],
                            score2=classification
                        ))
                    else:
                        pass
                else:
                    var_list.append(lookup)
                    class_list.append(row['Assigned Pathogenicity Score'])

        errors.write("Total variants reported more than once: {repeats}.\n".format(repeats=repeats))
        errors.write("Total conflicts: {conflicts}.\n".format(conflicts=conflicts))

