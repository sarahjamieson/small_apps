import re
from suds.client import Client
import sys
import socket
import urllib2


def get_gene(search_text):
    """Uses regular expressions to search for and extract "BRCA1" and "BRCA2" in a piece of text.

    Note: outputs an empty string if no gene is found.

    :param search_text: any string.
    :return: "BRCA1" or "BRCA2".
    """
    gene_search = re.search("BRCA1|BRCA2|CHD1|PTEN|STK11|TP53|ATM|BRIP1|CHEK2|PALB2|RAD51C|RAD51D", search_text)
    if gene_search:
        gene = gene_search.group(0)
    else:
        gene = ""

    return gene


def get_class(summary):
    """Takes the summary line of the Shire report and determines the variant classification based on keywords or phrases
    in the text.

    Uses regular expressions to look for all words/phrases associated with each classification type; these currently
    cover all variations found so far. Class 3 is search for twice as one summary found reported a class 3 as a
    "missense mutation" which would otherwise be classified as "U". If this was added to the original class 3 keywords,
    results reported as "pathogenic missense" would be incorrectly classified as a "3".

    Note:- Class 1 and 2 variants are not usually included in the report unless it is a variant update/review.

    :param summary: summary line of report (e.g. No evidence of a pathogenic mutation...)
    :return: the classification (1, 2, 3, 4, 5) or N for negative or U for unknown.
    """

    if re.search(r"unlikely", summary.lower()) > 0:
        var_class = "2"
    elif re.search(r"not\sclinically\simportant|benign|polymorphism", summary.lower()) > 0:
        var_class = "1"
    elif re.search(r"no\sevidence|no\sapparent\sevidence|normal|neg|-ve|no\smut|no\sother\svariants", summary.lower()) > 0:
        var_class = "N"
    elif re.search(r"likely\spathogenic|consistent", summary.lower()) > 0:
        var_class = "4"
    elif re.search(r"pathogenic|out-of-frame", summary.lower()) > 0:
        var_class = "5"
    elif re.search(r"uv|uncertain|missense\svariant|unclassified|unknown|variant|in-frame|heterozygous\s(deletion|duplication)", summary.lower()) > 0:
        var_class = "3"
    elif re.search(r"pathogenic|confirm|frameshift|nonsense|splice\ssite\smutation|deletion|splicesite|mutation",
                   summary.lower()) > 0:
        var_class = "5"
    elif re.search(r"missense", summary.lower()) > 0:
        var_class = "3"
    else:
        var_class = "U"

    return var_class


def get_dna_hgvs(summary):
    """Searches a given string and returns the first cDNA nomenclature within the string.

    The search string would typically be the summary line of a Shire report. If a variant has been found, the cDNA
    change will always be included in the summary line, unlike protein changes.

    Note: only the first pattern will be returned, so if there is >1 in the string these should be split beforehand.

    For variants with a classification other than 3, 4 and 5 and for variants where the HGVS cannot be determined,
    the cDNA HGVS will be a blank string.

    :param summary: search text e.g. summary line of report (No evidence of a pathogenic mutation...)
    :return: cDNA HGVS (dna_change).
    """

    c_search = re.search(r"c\.[^ ]+\b", summary)
    if c_search:
        dna_change = c_search.group(0)
    else:
        dna_change = ""
    return dna_change


def format_p_hgvs(protein_change):
    if "p." in protein_change:
        if ")" in protein_change and "(" not in protein_change:
            protein_change = protein_change[:2] + "(" + protein_change[2:]
        elif "(" in protein_change and ")" not in protein_change:
            protein_change += ")"
        elif "(" not in protein_change and ")" not in protein_change:
            protein_change = protein_change[:2] + "(" + protein_change[2:] + ")"
        elif "))" in protein_change:
            protein_change = protein_change[:-1]
        else:
            pass
        if "X" in protein_change:
            protein_change = protein_change.replace("X", "*")
    else:
        pass

    return protein_change


def get_protein_hgvs(summary, report):
    """Searches two strings and returns the first protein nomenclature within both strings.

    The search strings would typically be the summary line and full report from Shire. Protein changes are not always
    mentioned in the summary line so function first checks summary line and, if no "p." pattern is identified, it
    checks the full report.

    :param summary: search text e.g. summary line of report (No evidence of a pathogenic mutation...)
    :param report: other search text e.g. full report section.
    :return: protein HGVS (protein_change)
    """
    p_search_summary = re.search(r"p\.[^.,\s]+", summary)
    p_search_report = re.search(r"p\.[^.,\s]+", report)

    if p_search_summary:
        protein_change = p_search_summary.group(0)
    elif p_search_report:
        protein_change = p_search_report.group(0)
    else:
        protein_change = ""
    protein_change = format_p_hgvs(protein_change)

    return protein_change


def get_mlpa_result(summary):
    all_exons = re.findall(r'(?<!BRCA)\d+', summary)
    if "duplication" in summary:
        dup_or_del = "Dup"
    elif "deletion" in summary:
        dup_or_del = "Del"
    else:
        dup_or_del = "U"
    if len(all_exons) == 1:
        exon_range = all_exons[0].replace(" ", "")
        mlpa_result = dup_or_del + " Exon(s) " + exon_range
    elif len(all_exons) == 2:
        exon_range = all_exons[0].replace(" ", "") + "-" + all_exons[1].replace(" ", "")
        mlpa_result = dup_or_del + " Exon(s) " + exon_range
    elif not all_exons:
        mlpa_result = dup_or_del + " Entire Gene"
    else:
        raise Exception("Cannot determine exon range.")

    return mlpa_result


def get_genome_hgvs(genome_build, ref_seq, cdna_hgvs):
    """Checks Mutalyzer database for coordinate information for a given variant.

    :param genome_build: e.g. hg19
    :param ref_seq: e.g. NM_007294.2 (version required)
    :param cdna_hgvs: e.g. c.123G>C
    :return:
    """
    query = "{refseq}:{cdna}".format(refseq=ref_seq, cdna=cdna_hgvs)
    url = "https://mutalyzer.nl/services/?wsdl"

    client = Client(url, cache=None)
    response = client.service.numberConversion(genome_build, query)
    find_hgvs = re.search("\".*\"", str(response))
    if find_hgvs:
        genome_hgvs = find_hgvs.group(0)[1:-1]
    else:
        genome_hgvs = ""
    print genome_hgvs

    return genome_hgvs

