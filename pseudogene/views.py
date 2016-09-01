from django.shortcuts import render
from search import get_pseudogene_ids, create_fasta
from create_html import create_html
import os


def search_form(request):
    error = False
    zero_pseudo = False
    if 'gene' in request.GET:
        gene = request.GET['gene'].upper()
        if not gene:
            error = True
        else:
            # add if gene.html already exists, show straight away
            no_of_pseudogenes, pseudogene_dict, gene_seq = get_pseudogene_ids(gene)
            if no_of_pseudogenes == 0:
                zero_pseudo = True
            else:
                infile = create_fasta(no_of_pseudogenes, pseudogene_dict, gene_seq, gene)
                html_file = create_html(infile, no_of_pseudogenes, gene)
                html_page = '/home/shjn/PycharmProjects/small_apps/pseudogene/templates/pseudogene/%s' % html_file
                os.system("mv /home/shjn/PycharmProjects/small_apps/%s "
                          "/home/shjn/PycharmProjects/small_apps/pseudogene/templates/pseudogene/" % html_file)
                return render(request, 'pseudogene/results.html', {'html_page': html_page})

    return render(request, 'pseudogene/index.html', {'error': error, 'zero': zero_pseudo})
