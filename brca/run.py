from subprocess import Popen, PIPE

p = Popen(["python", "/home/shjn/PycharmProjects/brca_challenge/brca.py", "-i", "/home/shjn/PycharmProjects/brca_challenge/inputs/brca_1314.csv", "-o", "output", "-e", "errors.txt", "-p", "/media/sf_sarah_share/BRCA_poly_list.xls"], stdout=PIPE, stderr=PIPE)
out, err = p.communicate()
print out
print err

