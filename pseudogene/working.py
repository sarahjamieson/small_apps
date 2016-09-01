import re
import xml.etree.ElementTree as ET
from urllib2 import urlopen
from Bio import AlignIO
from biomart import BiomartServer
import json
import os

import requests
import sys

'''
def get_pseudogene_ids(gene):
    ensg = ""
    request_url = urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]+AND+"
        "Homo+sapiens[Organism]" % gene
    )
    tree = ET.parse(request_url)
    root = tree.getroot()
    pseudogene_ids = []
    pseudo_seqs = []
    pseudo_names = []

    for id_list in root.findall('IdList'):
        gene_id = id_list.find('Id').text
        gene_url = urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=related_functional_gene_%s[group]"
            % gene_id
        )
        pseudo_tree = ET.parse(gene_url)
        for parent in pseudo_tree.getiterator('IdList'):
            for child in parent:
                id_no = child.text
                pseudogene_ids.append(id_no)
        ensg_url = urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s&rettype=xml&retmode=text" % gene_id
        )
        ensg_tree = ET.parse(ensg_url)
        ensg_root = ensg_tree.getroot()
        for item in ensg_root.findall(
                'Entrezgene/Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_tag/Object-id/Object-id_str'
        ):
            if re.match("ENSG(.*)", item.text):
                ensg = item.text

    no_of_pseudogenes = len(pseudogene_ids)
    gene_seq = get_sequence(ensg)
    for pseudogene in pseudogene_ids:
        name = get_pseudogene_name(pseudogene)
        pseudo_names.append(name)
        start, end, chrom = get_gi(pseudogene)
        pseudo_seq = get_pseudo_seq(chrom, start, end)
        pseudo_seqs.append(pseudo_seq)

    pseudogene_dict = dict(zip(pseudogene_ids, zip(pseudo_names, pseudo_seqs)))
    print pseudogene_dict
    return no_of_pseudogenes, pseudogene_dict, gene_seq


def get_sequence(ensg):
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/%s?" % ensg
    r = requests.get(server+ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    gene_seq = r.text

    return gene_seq


def get_pseudogene_seq(pseudogene):
    pseudo_ensg = None
    seq = None
    server = "https://rest.ensembl.org"
    ext = "/xrefs/symbol/homo_sapiens/%s?" % pseudogene
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    for item in decoded:
        for key, value in item.iteritems():
            if key == "id" and re.match("ENSG(.*)", value):
                pseudo_ensg = value

    if pseudo_ensg is not None:
        seq = get_sequence(pseudo_ensg)

    return seq
'''

'''
def get_pseudogene_name(pseudo_id):
    pseudogene_name = None
    request_url = urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s" % pseudo_id)
    tree = ET.parse(request_url)
    root = tree.getroot()
    for item in root.findall('DocumentSummarySet/DocumentSummary/Name'):
        pseudogene_name = item.text
    return pseudogene_name


def create_fasta(pseudogene_dict, gene_seq, gene):
    with open("%s.txt" % gene, "w+") as fasta_file:
        for key, value in pseudogene_dict.iteritems():
            fasta_file.write(">%s\n%s\n" % (key, value))
        fasta_file.write(">%s\n%s\n" % (gene, gene_seq))

#no_of_pseudogenes, pseudogene_dict, gene_seq = get_pseudogene_ids("SDHD")
#create_fasta(pseudogene_dict, gene_seq, "SDHD")


def get_gi(pseudogene_id):
    seq_from = None
    seq_to = None
    chrom = None
    request_url = urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s&rettype=xml&retmode=text"
        % pseudogene_id
    )
    tree = ET.parse(request_url)
    root = tree.getroot()
    for item in root.find(
            'Entrezgene/Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval'
    ):
        if item.tag == "Seq-interval_from":
            seq_from = item.text
        elif item.tag == "Seq-interval_to":
            seq_to = item.text
    if seq_from < seq_to:
        start = seq_from
        end = seq_to
    else:
        start = seq_to
        end = seq_from
    for item in root.find(
        'Entrezgene/Entrezgene_source/BioSource/BioSource_subtype/SubSource'
    ):
        if item.tag == "SubSource_name":
            chrom = item.text
    return start, end, chrom


def get_pseudo_seq(chrom, start, end):
    pseudo_seq = None
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/%s:%s..%s?" % (chrom, start, end)
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    for key, value in decoded.iteritems():
        if key == "seq":
            pseudo_seq = value

    return pseudo_seq

# os.system("clustalo -i SDHD.txt -o SDHD_alignment.txt --outfmt=fa --force --resno")

alignment = AlignIO.read("SDHD_alignment.txt", "fasta")
total_length = len(alignment[0].seq)
x = 0
no_of_pseudogenes = 7
with open("templates/pseudogene/SDHD_test.html", "w+") as test_file:
    test_file.write("<!DOCTYPE html>\n<html lang=\"en\">\n<link rel=\"stylesheet\" "
                    "href=\"http://www.w3schools.com/lib/w3.css\">\n<head>\n<meta charset=\"UTF-8\">\n"
                    "<title>Results</title>\n<style>\nmark {background-color: yellow;}\n#align_table {display: block;}"
                    "\n#align_table td {display: inline-block;}\n</style>\n</head>\n<body>\n<h2>Alignment result</h2>\n"
                    "<p>There are %s pseudogenes in this gene.</p>\n<table id=\"align_table\">\n<tr>\n"
                    % no_of_pseudogenes)

    while x < total_length:
        base_alignment = alignment[:, x]
        ref_base = base_alignment[0]
        test_file.write("<td>\n<pre><font color=\"red\">%s</font>\n" % ref_base)
        for base in base_alignment[1:(no_of_pseudogenes + 1)]:
            if base == ref_base:
                test_file.write("%s\n" % base)
            elif base == "-":
                test_file.write("%s\n" % base)
            else:
                test_file.write("<mark>%s</mark>\n" % base)
        test_file.write("</pre>\n</td>\n")
        x += 1

    test_file.write("</tr>\n</table>\n</body>\n</html>")
'''
'''
alignment = AlignIO.read("SDHD_alignment.txt", "fasta")
total_length = len(alignment[0].seq)
x = 0
y = 0
no_of_pseudogenes = 7
with open("templates/pseudogene/SDHD_test.html", "w+") as test_file:
    test_file.write("<!DOCTYPE html>\n<html lang=\"en\">\n<link rel=\"stylesheet\" "
                    "href=\"http://www.w3schools.com/lib/w3.css\">\n<head>\n<meta charset=\"UTF-8\">\n"
                    "<title>Results</title>\n<style>\npre.space {line-height: 200%%}\nmark {background-color: yellow;}"
                    "\n#align_table {display: block;}\n#align_table td {display: inline-block;}\n</style>\n</head>\n"
                    "<body>\n<h2>Alignment result</h2>\n<p>There are %s pseudogenes in this gene.</p>\n"
                    % no_of_pseudogenes)
    while y < 555:
        id_value = 0
        sixty_bases = alignment[:, :60]
        bases = len(sixty_bases[0].seq)
        test_file.write("<table id=\"align_table\">\n<tr>\n<td>\n<pre id=\"space\">\n")
        while id_value < no_of_pseudogenes:
            test_file.write("%s\n" % sixty_bases[id_value].id)
            id_value += 1
        test_file.write("</pre>\n</td>\n")
        while x < bases:
            base_alignment = sixty_bases[:, x]
            ref_base = base_alignment[0]
            test_file.write("<td>\n<pre><font color=\"red\">%s</font>\n" % ref_base)
            for base in base_alignment[1:(no_of_pseudogenes + 1)]:
                if base == ref_base:
                    test_file.write("%s\n" % base)
                elif base == "-":
                    test_file.write("%s\n" % base)
                else:
                    test_file.write("<mark>%s</mark>\n" % base)
            test_file.write("</pre>\n</td>\n")
            x += 1
        test_file.write("</tr>\n</table>\n")
        # end table
        alignment = alignment[:, 60:]
        x = 0
        y += 1


# then make a table for each base_alignment, add column with names, then for base in base_alignment columns, then count

alignment = AlignIO.read("SDHD_alignment.txt", "fasta")
total_length = len(alignment[0].seq)
no_of_pseudogenes = 7
y = 0
a = 0
counts_dict = {}
while a < 7:
    counts_dict["%s" % a] = 0
    a += 1

with open("templates/pseudogene/SDHD_test.html", "w+") as test_file:
    test_file.write("<!DOCTYPE html>\n<html lang=\"en\">\n<link rel=\"stylesheet\" "
                    "href=\"http://www.w3schools.com/lib/w3.css\">\n<head>\n<meta charset=\"UTF-8\">\n"
                    "<title>Results</title>\n<style>\nmark {background-color: yellow;}"
                    "\n.table {display: table;  font-family: \"Courier New\"; font-size: medium; }\n.row { display: table-row; }\n.cell1 "
                    "{ display: table-cell; width: 200px; }\n.cell2 { display: table-cell; width: 1450px; }\n"
                    ".cell3 { display: table-cell; width: 200px; }\n</style>\n</head>\n<body>\n<h2>Alignment result"
                    "</h2>\n<p>There are %s pseudogenes in this gene.</p>\n" % no_of_pseudogenes)
    while y < 222:
        test_file.write("<div class=\"table\">\n")
        x = 0
        sixty_bases = alignment[:, :150]
        while x < no_of_pseudogenes:
            sequence = sixty_bases[x].seq
            count_of_bases = sum(sequence.count(c) for c in ("A", "T", "C", "G"))
            counts_dict["%s" % x] += count_of_bases
            test_file.write(
                "<div class=\"row\">\n<div class=\"cell1\">%s</div>\n<div class=\"cell2\">" % sixty_bases[x].id
            )
            for base in sequence:
                test_file.write("%s" % base)
            test_file.write("</font></div>\n<div class=\"cell3\">%s</div>\n</div>\n" % counts_dict["%s" % x])
            x += 1
        test_file.write("</div>\n<p><br/></p>")
        alignment = alignment[:, 150:]
        y += 1
'''

'''
def get_pseudogene_ids(gene):
    ensg = ""
    request_url = urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene Name]+AND+"
        "Homo+sapiens[Organism]" % gene
    )
    tree = ET.parse(request_url)
    root = tree.getroot()
    pseudogene_ids = []
    pseudo_seqs = []
    pseudo_names = []

    for id_list in root.findall('IdList'):
        print id_list.text
        gene_id = id_list.find('Id').text
        print gene_id

get_pseudogene_ids("BRCA1")
'''


def get_pseudogene_ids(gene):
    ensg = ""
    request_url = urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=%s[Gene%%20Name]+AND+Homo+sapiens[Organism]" % gene
    )
    tree = ET.parse(request_url)
    root = tree.getroot()
    pseudogene_ids = []
    pseudo_seqs = []
    pseudo_names = []
    for id_list in root.findall('IdList'):
        gene_id = id_list.find('Id').text
        gene_url = urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=related_functional_gene_%s[group]"
            % gene_id
        )
        pseudo_tree = ET.parse(gene_url)
        for parent in pseudo_tree.getiterator('IdList'):
            for child in parent:
                id_no = child.text
                pseudogene_ids.append(id_no)
        ensg_url = urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s&rettype=xml&retmode=text" % gene_id
        )
        ensg_tree = ET.parse(ensg_url)
        ensg_root = ensg_tree.getroot()
        for item in ensg_root.findall(
                'Entrezgene/Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag/Dbtag_tag/Object-id/Object-id_str'
        ):
            if re.match("ENSG(.*)", item.text):
                ensg = item.text

    no_of_pseudogenes = len(pseudogene_ids)
    pseudogene_dict = dict(zip(pseudogene_ids, zip(pseudo_names, pseudo_seqs)))
    gene_seq = get_sequence(ensg)
    for pseudogene in pseudogene_ids:
        name = get_pseudogene_name(pseudogene)
        pseudo_names.append(name)
        start, end, chrom = get_gi(pseudogene)
        pseudo_seq = get_pseudo_seq(chrom, start, end)
        pseudo_seqs.append(pseudo_seq)

    pseudogene_dict = dict(zip(pseudogene_ids, zip(pseudo_names, pseudo_seqs)))

    return no_of_pseudogenes, pseudogene_dict, gene_seq


def get_sequence(ensg):
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/%s?" % ensg
    r = requests.get(server+ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    gene_seq = r.text

    return gene_seq


def get_pseudogene_name(pseudo_id):
    pseudogene_name = None
    request_url = urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s" % pseudo_id)
    tree = ET.parse(request_url)
    root = tree.getroot()
    for item in root.findall('DocumentSummarySet/DocumentSummary/Name'):
        pseudogene_name = item.text
    return pseudogene_name


def get_gi(pseudogene_id):
    seq_from = None
    seq_to = None
    chrom = None
    request_url = urlopen(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%s&rettype=xml&retmode=text"
        % pseudogene_id
    )
    tree = ET.parse(request_url)
    root = tree.getroot()
    for item in root.find(
            'Entrezgene/Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval'
    ):
        if item.tag == "Seq-interval_from":
            seq_from = item.text
        elif item.tag == "Seq-interval_to":
            seq_to = item.text
    if seq_from < seq_to:
        start = seq_from
        end = seq_to
    else:
        start = seq_to
        end = seq_from
    for item in root.find(
        'Entrezgene/Entrezgene_source/BioSource/BioSource_subtype/SubSource'
    ):
        if item.tag == "SubSource_name":
            chrom = item.text
    return start, end, chrom


def get_pseudo_seq(chrom, start, end):
    pseudo_seq = None
    server = "https://rest.ensembl.org"
    ext = "/sequence/region/human/%s:%s..%s?" % (chrom, start, end)
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    for key, value in decoded.iteritems():
        if key == "seq":
            pseudo_seq = value

    return pseudo_seq


def create_fasta(pseudogene_dict, gene_seq, gene, no_of_pseudogenes):
    with open("%s.txt" % gene, "w+") as fasta_file:
        fasta_file.write(">%s_gene\n%s\n" % (gene, gene_seq))
        for key, value in pseudogene_dict.iteritems():
            fasta_file.write(">%s_%s\n%s\n" % (key, value[0], value[1]))
    if no_of_pseudogenes != 0:
        os.system("clustalo -i %s.txt -o %s_alignment.txt --outfmt=fa --force --resno" % (gene, gene))
    filename = "%s_alignment.txt" % gene

    return filename

no_of_pseudogenes, pseudogene_dict, gene_seq = get_pseudogene_ids("SOX9")
create_fasta(pseudogene_dict, gene_seq, "SOX9", no_of_pseudogenes)

