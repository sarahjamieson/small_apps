import re
import os


def check_for_known_errors(dataframe, error_file, indication):
    """Checks for and highlights errors previously identified in test data which can only be remedied manually.

    Possible errors:
    - Protein nomenclature that doesn't end with a closed bracket (indicates space in initial p. nomenclature)
    - Missing protein nomenclature or misleading/insufficient (often the case with synonymous variants)

    :param indication:
    :param error_file:
    :param dataframe: to check.
    :return:
    """
    cwd = os.getcwd()
    with open(error_file, "ab+") as errors:
        dataframe = dataframe.reset_index(drop=True)
        for row_index, row in dataframe.iterrows():
            # If spaces in p. nomenclature then will not have ended with a ")".
            if row["Protein Impact"] and not row["Protein Impact"].endswith(")"):
                if "NHS Number" in dataframe.columns:
                    errors.write("{nhs_no}\tError in protein nomenclature.\t{indication}\n".format(
                        nhs_no=row["NHS Number"], indication=indication)
                    )
                else:
                    errors.write("{cdna}\tError in protein nomenclature.\t{indication}\n".format(
                        cdna=row[3], indication=indication)
                    )
            '''
            if row["cDNA Change"] and not re.match("Del|Dup", row["cDNA Change"]) and not row["Genomic Change"]:
                if "NHS Number" in dataframe.columns:
                    errors.write("{nhs_no}\tMissing genomic HGVS nomenclature. Check validity of cDNA HGVS."
                                 "\t{indication}\n".format(nhs_no=row["NHS Number"], indication=indication))
                else:
                    errors.write("{cdna}\tMissing genomic HGVS nomenclature. Check validity of cDNA HGVS.\t"
                                 "{indication}\n".format(cdna=row[3], indication=indication))
            '''
            # Checks for missing or misleading/insufficient p. nomenclature. Splicing and MLPA results excluded.
            no_protein_hgvs = ["+", "-", "DEL", "DUP"]
            no_protein_variants = 0
            for item in no_protein_hgvs:
                if item in row["cDNA Change"].upper():
                    no_protein_variants += 1
            if row["cDNA Change"]:
                if re.match("\.\D", row["cDNA Change"]):
                    if "NHS Number" in dataframe.columns:
                        errors.write("{nhs_no}\tUnexpected letter(s) in cDNA nomenclature.\t{indication}\n".format(
                            nhs_no=row["NHS Number"], indication=indication))
                    else:
                        errors.write("{cdna}\tUnexpected letter(s) in cDNA nomenclature.\t{indication}\n".format(
                            cdna=row[3], indication=indication)
                        )
                if no_protein_variants == 0 and not row["Protein Impact"]:
                    if "NHS Number" in dataframe.columns:
                        errors.write("{nhs_no}\tMissing protein nomenclature.\t{indication}\n".format(
                            nhs_no=row["NHS Number"], indication=indication)
                        )
                    else:
                        errors.write("{cdna}\tMissing protein nomenclature.\t{indication}\n".format(
                            cdna=row[3], indication=indication)
                        )
                elif re.match("p.(\?|\(\?\)|\(=\))", row["Protein Impact"]) or ";" in row["Protein Impact"]:
                    if "NHS Number" in dataframe.columns:
                        errors.write("{nhs_no}\tIf this is a silent change, this should be in the format e.g. p."
                                     "(Leu54=).\t{indication}\n".format(nhs_no=row["NHS Number"], indication=indication))
                    else:
                        errors.write("{cdna}\tIf this is a silent change, this should be in the format e.g. p.(Leu54=)."
                                     "\t{indication}\n".format(cdna=row[3], indication=indication))
            else:
                if "NHS Number" in dataframe.columns:
                    errors.write("{nhs_no}\tNo cDNA Change detected.\t{indication}\n".format(
                        nhs_no=row["NHS Number"], indication=indication)
                    )
                else:
                    errors.write("{g}\tNo cDNA Change detected.\t{indication}\n".format(
                        g=row["Genomic Change"], indication=indication)
                    )
            if row["Gene"] is None:
                if "NHS Number" in dataframe.columns:
                    errors.write("{nhs_no}\tNo gene name detected.\t{indication}\n".format(
                        nhs_no=row["NHS Number"], indication=indication)
                    )
                else:
                    errors.write("{cdna}\tNo gene name detected.\t{indication}\n".format(
                        cdna=row["cDNA Change"], indication=indication)
                    )

            # no space after c. so may include p. too. Add full stop followed by letter.
            if "," in row["cDNA Change"]:
                cdna_split = row["cDNA Change"].split(",")
                cdna = cdna_split[0]
                dataframe = dataframe.set_value(row_index, "cDNA Change", cdna)
    return dataframe
