import re
import sys
import csv
import datetime
import pandas as pd
from polys import get_brca_polys
from check_conflicts import check_conflicts
from known_errors import check_for_known_errors
from further_studies import CheckAgainstFurtherStudies
from get_attributes import get_gene, get_class, get_protein_hgvs, get_dna_hgvs, get_mlpa_result, format_p_hgvs, \
    get_genome_hgvs
import os


def run_brca_extraction(input_csv, output_prefix, error_file, poly_file):
    # Declarations
    cwd = os.getcwd()
    all_classes_in_row = []
    poly_classes = ["1", "2"]
    var_classes = ["3", "4", "5"]
    full_col_list = ["NHS Number", "Forename(s)", "Surname", "DOB", "Postcode", "Molecular Testing Type", "Gene",
                     "RefSeq Transcript ID", "Genomic Change", "cDNA Change", "Protein Impact",
                     "Assigned Pathogenicity Score", "Sample Report Date"]
    reduced_col_list = ["Gene", "RefSeq Transcript ID", "Genomic Change", "cDNA Change", "Protein Impact"]
    neg_output_df = pd.DataFrame(columns=full_col_list)
    all_results_df = pd.DataFrame(columns=full_col_list)
    main_output_df = pd.DataFrame(columns=full_col_list)
    poly_output_df = pd.DataFrame(columns=reduced_col_list)
    unknown_output_df = pd.DataFrame(columns=full_col_list)

########################################################################################################################

    # Extract variant classifications determined by RNA studies and variant updates into two dictionaries.
    print "Collating all RNA studies and variant update/review data...",
    sys.stdout.flush()

    check_others = CheckAgainstFurtherStudies(input_csv)
    rna_studies = check_others.check_rna_studies()
    updates = check_others.check_variant_updates()
    print "[DONE]"

########################################################################################################################

    print "Extracting variants...",
    sys.stdout.flush()

    with open(input_csv, "rb") as infile:
        reader = csv.reader(infile, delimiter='\t')
        next(reader, None)  # skips header row
        for row in reader:

            all_classes_in_row = []
            summary = row[11]
            report = row[13]

            if row[0] and row[14] != "PP" and not re.match(
                    "Tumour\sstudies|RNA|Variant\supdate|Variant\sreview", row[8]
            ) and not re.match("VARIANT\sREVIEW", summary):
                with open(error_file, "ab+") as errors:
                    if "heterozygous missense mutation" in summary.lower():
                        errors.write(
                            "{nhs_no}\tPotential incorrect classification.\tSee \"main\" file.\n".format(nhs_no=row[0]))
                    if "amended" in summary.lower():
                        errors.write(
                            "{nhs_no}\tHas had amended report, please check for duplicates or nomenclature changes.\t"
                            "Check results.\n".format(nhs_no=row[0]))
                    if "but" in summary.lower():
                        errors.write("{nhs_no}\tPotential missed variant.\tSee \"main\" file.\n".format(nhs_no=row[0]))
                    if re.match("(alt|disrupt|ab)(.*)(splicing)", summary):
                        errors.write("{nhs_no}\tUnspecified RNA studies.\tCheck results.\n".format(nhs_no=row[0]))
                    if "loss of heterozygosity" in summary.lower():
                        errors.write("{nhs_no}\tUnspecified tumour studies.\tCheck results.\n".format(nhs_no=row[0]))

                # Calculate number of summary "segments" and number of variants.
                summaries = re.split(r"\.\s|\.\\n", summary)
                variants = re.findall(r"c\.[^ ]+\b", summary)
                mlpa_variants = re.findall(r"(pathogenic|heterozygous|out-of-frame)\s(deletion|duplication)", summary.lower())

                # (i) Report with single MLPA result.
                if len(summaries) == 1 and len(mlpa_variants) == 1 and len(variants) == 0:
                    var_class = get_class(summary)

                    if var_class == "N":
                        var_df = pd.DataFrame([[row[0], row[3], row[4], row[2], row[15], row[8], "", "", "", "", "", "N", row[9]]],
                                              columns=full_col_list)
                        all_results_df = all_results_df.append(var_df)
                        all_classes_in_row.append(var_class)
                    else:
                        if var_class == "3":
                            if "pathogenic deletion" in report.lower() or "pathogenic duplication" in report.lower():
                                var_class = "5"

                        c_hgvs = get_mlpa_result(summary)
                        gene = get_gene(summary)

                        if gene == "BRCA1" or gene == "BRCA2":
                            if gene == "BRCA1":
                                refseq = "NM_007294.3"
                            elif gene =="BRCA2":
                                refseq = "NM_000059.3"

                            var_df = pd.DataFrame([[row[0], row[3], row[4], row[2], row[15], row[8], gene, refseq, "",
                                                    c_hgvs, "", var_class, row[9]]], columns=full_col_list)
                            all_results_df = all_results_df.append(var_df)
                            all_classes_in_row.append(var_class)

                # (ii) MLPA and sequencing (c.) result.
                elif len(mlpa_variants) == 1 and len(variants) == 1:

                    if len(summaries) > 1:
                        summary_split = summaries
                    else:
                        summary_split = summary.split("and")

                    for item in summary_split:
                        if "c." in item:
                            var_class = get_class(item)
                            c_hgvs = variants[0]
                            p_hgvs = get_protein_hgvs(item, report)
                            gene = get_gene(item)

                        else:
                            var_class = get_class(item)
                            c_hgvs = get_mlpa_result(item)
                            p_hgvs = ""
                            gene = get_gene(item)

                        if gene == "BRCA1" or gene == "BRCA2":
                            if gene == "BRCA1":
                                refseq = "NM_007294.3"
                            elif gene == "BRCA2":
                                refseq = "NM_000059.3"
                            g_hgvs = get_genome_hgvs("hg19", refseq, c_hgvs)
                            var_df = pd.DataFrame([[row[0], row[3], row[4], row[2], row[15], row[8], gene, refseq, g_hgvs,
                                                    c_hgvs, p_hgvs, var_class, row[9]]], columns=full_col_list)
                            all_results_df = all_results_df.append(var_df)
                            all_classes_in_row.append(var_class)

                # (iii) Sequencing result, one conclusion.
                elif len(summaries) == 1 and len(variants) < 2:
                    var_class = get_class(summary)

                    if var_class != "N" and var_class != "U":
                        c_hgvs = get_dna_hgvs(summary)
                        p_hgvs = get_protein_hgvs(summary, report)
                        gene = get_gene(summary)
                    else:
                        gene = ""
                        c_hgvs = ""
                        p_hgvs = ""

                    if var_class == "N":
                        var_df = pd.DataFrame(
                            [[row[0], row[3], row[4], row[2], row[15], row[8], gene, "", "", c_hgvs, p_hgvs, var_class, row[9]]],
                            columns=full_col_list
                        )
                        all_results_df = all_results_df.append(var_df)
                        all_classes_in_row.append(var_class)
                    elif gene == "BRCA1" or gene == "BRCA2":
                        if gene == "BRCA1":
                            refseq = "NM_007294.3"
                        elif gene == "BRCA2":
                            refseq = "NM_000059.3"

                        g_hgvs = get_genome_hgvs("hg19", refseq, c_hgvs)
                        var_df = pd.DataFrame(
                            [[row[0], row[3], row[4], row[2], row[15], row[8], gene, refseq, g_hgvs, c_hgvs, p_hgvs,
                              var_class, row[9]]], columns=full_col_list
                        )
                        all_results_df = all_results_df.append(var_df)
                        all_classes_in_row.append(var_class)

                # (iv) Two sequencing results in one sentence.
                elif len(summaries) == 1 and len(variants) > 1:
                    summary_split = summary.split('and')
                    var_class1 = get_class(summary_split[0])
                    var_class2 = get_class(summary_split[1])

                    if var_class2 == "U":
                        var_class2 = var_class1
                    elif var_class1 == "N":
                        var_class2 = var_class1
                    elif var_class2 == "N":
                        var_class1 = var_class2

                    if var_class1 == "N" and var_class2 == "N":
                        c_hgvs = ""
                        p_hgvs = ""
                        gene = ""
                        g_hgvs = ""
                        refseq = ""
                        var_df = pd.DataFrame(
                            [[row[0], row[3], row[4], row[2], row[15], row[8], gene, refseq, g_hgvs, c_hgvs, p_hgvs,
                              var_class, row[9]]], columns=full_col_list
                        )
                        all_results_df = all_results_df.append(var_df)
                        all_classes_in_row.append(var_class)
                    else:
                        c_hgvs_1 = variants[0]
                        c_hgvs_2 = variants[1]
                        p_search_summary = re.findall(r"p\.[^.,\s]+", summary)
                        p_search_report = re.findall(r"p\.[^.,\s]+", report)

                        if p_search_summary:
                            p_search = p_search_summary
                        elif p_search_report:
                            p_search = p_search_report
                        else:
                            raise Exception("No \"p.\" pattern could be found in the text.")

                        if "-" in c_hgvs_1 or "+" in c_hgvs_1:
                            p_hgvs_1 = ""
                            p_hgvs_2 = p_search[0]
                        elif "-" in c_hgvs_2 or "+" in c_hgvs_2:
                            p_hgvs_1 = p_search[0]
                            p_hgvs_2 = ""
                        else:
                            p_hgvs_1 = p_search[0]
                            if len(p_search) > 1:
                                p_hgvs_2 = p_search[1]
                            else:
                                p_hgvs_2 = ""

                        gene1 = get_gene(summary_split[0])
                        gene2 = get_gene(summary_split[1])

                        if not gene1:
                            gene1 = gene2
                        elif not gene2:
                            gene2 = gene1
                        else:
                            pass

                        p_hgvs_1 = format_p_hgvs(p_hgvs_1)
                        p_hgvs_2 = format_p_hgvs(p_hgvs_2)

                        if gene1 == "BRCA1" or gene1 == "BRCA2":
                            if gene1 == "BRCA1":
                                refseq = "NM_007294.3"
                            elif gene1 == "BRCA2":
                                refseq = "NM_000059.3"
                            g_hgvs_1 = get_genome_hgvs("hg19", refseq, c_hgvs_1)
                            var_df = pd.DataFrame(
                                [[row[0], row[3], row[4], row[2], row[15], row[8], gene1, refseq, g_hgvs_1, c_hgvs_1,
                                  p_hgvs_1, var_class1, row[9]]], columns=full_col_list
                            )
                            all_results_df = all_results_df.append(var_df)
                            all_classes_in_row.append(var_class)

                        if gene2 == "BRCA1" or gene2 == "BRCA2":
                            if gene2 == "BRCA1":
                                refseq = "NM_007294.3"
                            elif gene2 == "BRCA2":
                                refseq = "NM_000059.3"

                            g_hgvs_2 = get_genome_hgvs("hg19", refseq, c_hgvs_2)
                            var2_df = pd.DataFrame(
                                [[row[0], row[3], row[4], row[2], row[15], row[8], gene2, refseq, g_hgvs_2, c_hgvs_2,
                                  p_hgvs_2, var_class2, row[9]]], columns=full_col_list
                            )
                            all_results_df = all_results_df.append(var2_df)
                            all_classes_in_row.append(var_class)

                # (v) Two sentences, reporting negative result.
                elif len(summaries) > 1 and len(variants) == 0:
                    var_df = pd.DataFrame([[row[0], row[3], row[4], row[2], row[15], row[8], "", "", "", "", "", "N", row[9]]],
                                          columns=full_col_list)
                    all_results_df = all_results_df.append(var_df)
                    all_classes_in_row.append(var_class)

                # (vi) Two sentences, one sequencing result reported.
                elif len(summaries) > 1 and len(variants) == 1:
                    summary_split = summary.split(". ")

                    for item in summary_split:
                        var_class = get_class(item)
                        if var_class != "U" and var_class != "N":
                            c_hgvs = variants[0]
                            p_hgvs = get_protein_hgvs(summary, report)
                            gene = get_gene(item)

                            if gene == "BRCA1" or gene == "BRCA2":
                                if gene == "BRCA1":
                                    refseq = "NM_007294.3"
                                elif gene == "BRCA2":
                                    refseq = "NM_000059.3"

                                g_hgvs = get_genome_hgvs("hg19", refseq, c_hgvs)
                                var_df = pd.DataFrame([[row[0], row[3], row[4], row[2], row[15], row[8], gene, refseq,
                                                        g_hgvs, c_hgvs, p_hgvs, var_class, row[9]]], columns=full_col_list)
                                all_results_df = all_results_df.append(var_df)
                                all_classes_in_row.append(var_class)

                # (vii) Two sequencing results in two sentences.
                elif len(summaries) == 2 and len(variants) == 2:
                    summary_split = re.split(r"\.\s|\\n", summary)
                    print summary_split
                    for item in summary_split:
                        gene = get_gene(item)
                        item_index = summary_split.index(item)
                        print item_index
                        if gene == "BRCA1" or gene == "BRCA2":
                            var_class = get_class(item)
                            print variants
                            c_hgvs = variants[item_index]
                            if "-" in c_hgvs or "+" in c_hgvs:
                                p_hgvs = ""
                            else:
                                p_hgvs = re.search(r"p\.[^.,\s]+", item)
                                if p_hgvs:
                                    p_hgvs = p_hgvs.group(0)
                                else:
                                    p_hgvs = ""
                                p_hgvs = format_p_hgvs(p_hgvs)

                            if gene == "BRCA1":
                                refseq = "NM_007294.3"
                            elif gene == "BRCA2":
                                refseq = "NM_000059.3"

                            g_hgvs = get_genome_hgvs("hg19", refseq, c_hgvs)
                            var_df = pd.DataFrame([[row[0], row[3], row[4], row[2], row[15], row[8], gene, refseq, g_hgvs,
                                                    c_hgvs, p_hgvs, var_class, row[9]]], columns=full_col_list)
                            all_results_df = all_results_df.append(var_df)
                            all_classes_in_row.append(var_class)

                else:
                    with open(error_file, "ab+") as errors:
                        errors.write(
                            "{nhs_no}\tUnable to extract results from summary line.\tSee original file.\n".format(nhs_no=row[0]))

                with open(error_file, "ab+") as errors:
                    if "N" in all_classes_in_row and (
                                        "3" in all_classes_in_row or "4" in all_classes_in_row or "5" in all_classes_in_row
                    ):
                        errors.write("{nhs_no}\tBoth a positive and negative result has been found for patient.\t"
                                     "See \"main\" file.\n".format(nhs_no=row[0]))

                # Look for Class 2 variants only mentioned in the report.
                if "For information" in report:
                    split_report = report.split("For information")
                    fio_section = split_report[1]
                    fio_section_split = fio_section.split(". ")
                    fio_list = []
                    for sentence in fio_section_split:
                        if "in cis" not in sentence and "in trans" not in sentence:
                            fio_list.append(sentence)
                    fio_section = ". ".join(fio_list)
                    dna_changes = re.findall(r"c\.[^ ]+\b", fio_section)
                    protein_changes = re.findall(r"p\.[^.,\s]+", fio_section)
                    all_genes = re.findall(
                        "BRCA1|BRCA2|CHD1|PTEN|STK11|TP53|ATM|BRIP1|CHEK2|PALB2|RAD51C|RAD51D", fio_section
                    )

                    for item in dna_changes:
                        index = dna_changes.index(item)
                        if len(dna_changes) == len(protein_changes):
                            protein = format_p_hgvs(protein_changes[index])
                        elif len(protein_changes) == 0:
                            protein = ""
                        else:
                            if "+" in item or "-" in item:
                                protein = ""
                            else:
                                protein = format_p_hgvs(protein_changes[0])
                        if "A" not in protein and "T" not in protein and "G" not in protein and "C" not in protein:
                            protein = ""

                        if len(dna_changes) == len(all_genes):
                            gene = all_genes[index]
                        elif len(all_genes) == 0:
                            gene = ""
                        else:
                            gene = all_genes[0]

                        if gene == "BRCA1":
                            refseq = "NM_007294.3"
                        elif gene == "BRCA2":
                            refseq = "NM_000059.3"
                        if gene == "BRCA1" or gene == "BRCA2":
                            g_hgvs = get_genome_hgvs("hg19", refseq, item)
                            var_df = pd.DataFrame(
                                [[row[0], row[3], row[4], row[2], row[15], row[8], gene, refseq, g_hgvs, item, protein,
                                  "2", row[9]]], columns=full_col_list
                            )
                            all_results_df = all_results_df.append(var_df)
    print "[DONE]"

    ########################################################################################################################

    # Checks if any variants are in the RNA studies or variant updates dictionaries.
    print "Checking variants against RNA studies and variant updates/reviews...",
    sys.stdout.flush()

    all_results_df = all_results_df.reset_index(drop=True)

    for row_index, row in all_results_df.iterrows():
        cdna = row["cDNA Change"]
        gene = row["Gene"]
        report_date = datetime.datetime.strptime(row["Sample Report Date"], "%d/%m/%Y")
        combined_lookup = "{cdna}&{gene}".format(cdna=cdna, gene=gene)

        if combined_lookup in rna_studies:
            rna_result = rna_studies.get(combined_lookup).split("&")
            rna_class = rna_result[0]
            rna_date = datetime.datetime.strptime(rna_result[1], "%d/%m/%Y")
            if rna_date > report_date:
                all_results_df = all_results_df.set_value(row_index, "Assigned Pathogenicity Score", rna_class)

        if combined_lookup in updates:
            update_result = updates.get(combined_lookup).split("&")
            new_class = update_result[0]
            update_date = datetime.datetime.strptime(update_result[1], "%d/%m/%Y")
            if update_date > report_date:
                all_results_df = all_results_df.set_value(row_index, "Assigned Pathogenicity Score", new_class)

    print "[DONE]"

    ########################################################################################################################

    # Extracts recorded polys and adds to final poly output.
    print "Getting additional polys...",
    sys.stdout.flush()

    poly_df = get_brca_polys(poly_file)
    poly_output_df = poly_output_df.append(poly_df)

    print "[DONE]"

    ########################################################################################################################

    # Splits all variants by class into different dataframes.
    all_results_df = all_results_df.drop_duplicates()
    for row_index, row in all_results_df.iterrows():
        if row[11] in var_classes:
            row_df = pd.DataFrame(
                [[row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12]]],
                columns=full_col_list
            )
            main_output_df = main_output_df.append(row_df)
        elif row[11] == "N":
            row_df = pd.DataFrame(
                [[row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12]]],
                columns=full_col_list
            )
            neg_output_df = neg_output_df.append(row_df)
        elif row[11] in poly_classes:
            row_df = pd.DataFrame(
                [[row[6], row[7], row[8], row[9], row[10]]],
                columns=reduced_col_list
            )
            poly_output_df = poly_output_df.append(row_df)
        else:
            row_df = pd.DataFrame(
                [[row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12]]],
                columns=full_col_list
            )
            unknown_output_df = unknown_output_df.append(row_df)

    ########################################################################################################################

    # Checks each dataframe for known, unfixable errors.
    print "Checking data for known errors...",
    sys.stdout.flush()

    main_output_df = check_for_known_errors(main_output_df, error_file, "See \"main\" file.")
    # neg_output_df = check_for_known_errors(neg_output_df, args.error_file, "See \"neg\" file.")
    poly_output_df = check_for_known_errors(poly_output_df, error_file, "See \"poly\" file.")

    print "[DONE]"

    ########################################################################################################################

    # Checks for classification conflicts.
    print "Checking for classification conflicts...",
    sys.stdout.flush()

    check_conflicts(all_results_df, error_file)

    print "[DONE]"

    ########################################################################################################################

    # Write all dataframes to Excel tabs (change to CSV at end).
    print "Writing output to file...",
    sys.stdout.flush()

    main_output_df = main_output_df.drop_duplicates()
    main_output_df = main_output_df.sort_values(["NHS Number"])
    main_output_df.to_csv("{cwd}/brca/extraction/outputs/{prefix}_main.csv".format(cwd=cwd, prefix=output_prefix), index=False)

    if not poly_output_df.empty:
        poly_output_df = poly_output_df.drop_duplicates()
        poly_output_df.to_csv("{cwd}/brca/extraction/outputs/{prefix}_polys.csv".format(cwd=cwd, prefix=output_prefix), index=False)

    neg_output_df.to_csv("{cwd}/brca/extraction/outputs/{prefix}_neg.csv".format(cwd=cwd, prefix=output_prefix), index=False)

    if not unknown_output_df.empty:
        unknown_output_df.to_csv("{cwd}/brca/extraction/outputs/{prefix}_unknowns.csv".format(cwd=cwd, prefix=output_prefix), index=False)

    print "[DONE]"

    print "Data extraction complete."
