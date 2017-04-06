import os
import re
import pandas as pd
from bs4 import BeautifulSoup
import sys
import requests
import urllib2, cookielib
import time

# input: variant
# add variant to file
# run polyphen on file
# get results
# output values


def create_polyphen_input(chrom, pos, ref, alt):
    input_file = "%s_%s.txt" % (chrom, pos)
    with open(input_file, "w+") as i:
        i.write("%s:%s %s/%s\n" % (chrom, pos, ref, alt))
    i.close()

    return input_file


def run_polyphen(input_file):
    os.system("curl "
              "-F _ggi_project=PPHWeb2 "
              "-F _ggi_origin=query "
              "-F _ggi_target_pipeline=1 "
              "-F MODELNAME=HumDiv "  # HUMDIV or HUMVAR
              "-F UCSCDB=hg19 "  # genome
              "-F SNPFUNC=m "  # functional SNP category, c or m
              "-F _ggi_batch_file=@%s "  # filename of variants, format "e.g. chr1:1158631 A/C/G/T"
              "-D - http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi "
              "-c %s.humdiv.cookies.txt "
              "> %s.humdiv.txt" % (input_file, input_file[:-4], input_file[:-4]))

    os.system("curl "
              "-F _ggi_project=PPHWeb2 "
              "-F _ggi_origin=query "
              "-F _ggi_target_pipeline=1 "
              "-F MODELNAME=HumVar "  # HUMDIV or HUMVAR
              "-F UCSCDB=hg19 "  # genome
              "-F SNPFUNC=m "  # functional SNP category, c or m
              "-F _ggi_batch_file=@%s "  # filename of variants, format "e.g. chr1:1158631 A/C"
              "-D - http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi "
              "-c %s.humvar.cookies.txt"
              "> %s.humvar.txt" % (input_file, input_file[:-4], input_file[:-4]))
    humdiv_file = "%s.humdiv.txt" % input_file[:-4]
    humvar_file = "%s.humvar.txt" % input_file[:-4]

    return humdiv_file, humvar_file


def get_polyphen_results(hum_file):
    session_id = None
    opened_file = file(hum_file, "r").read()
    for word in opened_file.split():
        if re.match("polyphenweb2=.{40,};", word):
            session_id = word[13:53]
        else:
            pass
    if session_id is not None:
        print session_id
        a = 0
        polyphen_dict = {}
        cj = cookielib.MozillaCookieJar('%s.cookies.txt' % hum_file[:-4])
        cj.load()
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
        opener.addheaders = [('User-Agent', 'Mozilla/5.0')]
        home = opener.open("http://genetics.bwh.harvard.edu/ggi/pph2/%s/1/pph2-short.txt" % session_id)
        s = home.read()
        with open("hello.txt", "w+") as f:
            f.write(s)
        f.close()
        polyphen_df = pd.read_table("hello.txt", sep='\t', index_col=False)
        polyphen_df = polyphen_df.rename(columns=lambda x: x.strip())
        polyphen_df = polyphen_df.dropna()

        for row_index, row in polyphen_df.iterrows():
            polyphen_dict[0] = {
                'up_acc': row['acc'].strip(),
                'aa1': row['aa1'].strip(),
                'pos': row['pos'],
                'aa2': row['aa2'].strip(),
                'pred': row['prediction'].strip(),
                'prob_score': row['pph2_prob'],
                'sens': row['pph2_TPR'],
                'fpr': row['pph2_FPR']
            }
            a += 1

        return polyphen_dict

    else:
        print "No session ID detected"

chrom = 'chr13'
pos = '28608222'
ref = 'A'
alt = 'T'
filename = create_polyphen_input(chrom, pos, ref, alt)
humdiv, humvar = run_polyphen(filename)
time.sleep(30)
div_dict = get_polyphen_results(humdiv)
print(div_dict)
var_dict = get_polyphen_results(humvar)
print(var_dict)

