import re
import csv
from get_attributes import get_gene, get_dna_hgvs, get_class


class CheckAgainstFurtherStudies(object):
    def __init__(self, csvfile):
        self.csv = csvfile

    def check_rna_studies(self):
        with open(self.csv, "rb") as infile:
            reader = csv.reader(infile, delimiter='\t')
            next(reader, None)
            rna_studies = {}
            for row in reader:
                reason = row[8]
                summary = row[11]
                if "RNA studies" in reason:
                    cdna = re.search(r"c\.[^ ]+\b", summary)
                    if cdna:
                        cdna = cdna.group(0)
                    else:
                        cdna = ""
                    if "No evidence" in summary or "not" in summary:
                        var_class = "3"
                    elif "causes skipping" in summary:
                        var_class = "5"
                    else:
                        var_class = "U"
                    gene = get_gene(summary)
                    rna_studies["{cdna}&{gene}".format(cdna=cdna, gene=gene)] = "{var}&{date}".format(
                        var=var_class, date=row[9]
                    )

        return rna_studies

    def check_variant_updates(self):
        with open(self.csv, "rb") as infile:
            reader = csv.reader(infile, delimiter='\t')
            next(reader, None)
            updates = {}
            for row in reader:
                reason = row[8]
                summary = row[11]
                if "Variant Review" in reason or "Variant update" in reason:
                    var_class = get_class(summary)
                    cdna = get_dna_hgvs(summary)
                    gene = get_gene(summary)
                    updates["{cdna}&{gene}".format(cdna=cdna, gene=gene)] = "{var}&{date}".format(
                        var=var_class, date=row[9]
                    )

        return updates
