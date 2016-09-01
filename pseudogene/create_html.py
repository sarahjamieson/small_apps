from Bio import AlignIO
import os


def create_html(infile, no_of_pseudogenes, gene):
    alignment = AlignIO.read(infile, "fasta")
    total_length = len(alignment[0].seq)
    tables_required = total_length / 150
    y = 0
    a = 0
    counts_dict = {}
    while a < (no_of_pseudogenes + 1):
        counts_dict["%s" % a] = 0
        a += 1

    with open("%s.html" % gene, "w+") as test_file:
        test_file.write("<!DOCTYPE html>\n<html lang=\"en\">\n<link rel=\"stylesheet\" "
                        "href=\"http://www.w3schools.com/lib/w3.css\">\n<head>\n<meta charset=\"UTF-8\">\n<title>"
                        "Results</title>\n<style>\nmark {background-color: yellow;}\n.table {display: table; "
                        "font-family: \"Courier New\"; font-size: medium; }\n.cell1 { display: table-cell; width: "
                        "200px; }\n.cell2 { display: table-cell; width: 1450px; }\n.cell3 { display: table-cell; "
                        "width: 200px; }\n</style>\n</head>\n<body>\n<h2>Alignment result</h2>\n<p>There are %s "
                        "pseudogenes in this gene.</p>\n" % no_of_pseudogenes)
        while y < tables_required:
            test_file.write("<div class=\"table\">\n")
            x = 0
            sixty_bases = alignment[:, :150]
            while x < (no_of_pseudogenes + 1):
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
    html_file = "%s.html" % gene

    os.system("rm %s.txt" % gene)
    os.system("rm %s_alignment.txt")

    return html_file
