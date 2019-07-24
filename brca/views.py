from django.shortcuts import render, HttpResponseRedirect, HttpResponse
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.core.urlresolvers import reverse_lazy
from brca.forms import UploadFileForm
import os
from brca.extraction.brca import run_brca_extraction
import StringIO
import zipfile
import time
import socket
import urllib2


# http://stackoverflow.com/questions/15084597/django-error-message-for-login-form - try this for login forms
def user_login(request):
    logout(request)
    if request.method == 'POST':
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(username=username, password=password)
        if user:
            if user.is_active:
                login(request, user)
                return HttpResponseRedirect(reverse_lazy('home'))
            else:
                messages.ERROR(request, "Disabled account.")
        else:
            messages.ERROR(request, "Invalid login details.")
    return render(request, 'brca/login.html')


def user_logout(request):
    logout(request)
    return HttpResponseRedirect(reverse_lazy('home'))


def handle_uploaded_file(filename):
    """Uploads file in chunks.
        :param filename: file to upload.
        :return filepath: location of the file which has been uploaded.
    """
    filepath = os.path.join('brca/extraction/inputs/', filename.name)
    with open(filepath, 'wb') as destination:  # wb = file opened for writing in binary mode.
        for chunk in filename.chunks():
            destination.write(chunk)
    return filepath


def download(request):
    filepaths = ["brca/extraction/outputs/brca_output_main.csv",
                 "brca/extraction/outputs/brca_output_polys.csv",
                 "brca/extraction/outputs/brca_errors.txt"]
    zip_subdir = "BrcaExtraction_%s" % time.strftime("%d%m%Y")
    zip_filename = "%s.zip" % zip_subdir
    s = StringIO.StringIO()
    zf = zipfile.ZipFile(s, "w")
    for fpath in filepaths:
        fdir, fname = os.path.split(fpath)
        zip_path = os.path.join(zip_subdir, fname)
        zf.write(fpath, zip_path)
    zf.close()
    response = HttpResponse(s.getvalue(), content_type="application/x-zip-compressed")
    response['Content-Disposition'] = 'attachment; filename=%s' % zip_filename
    return response


@login_required(login_url=reverse_lazy('login'))
def home(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            filepath = handle_uploaded_file(request.FILES['file'])
            try:
                run_brca_extraction(filepath, "brca_output", "errors.txt", "/media/sf_sarah_share/BRCA_poly_list.xls")
                return render(request, 'brca/home.html', {'form': form, 'complete': True})
            except (socket.error, urllib2.URLError):
                return render(request, 'brca/home.html', {'form': form, 'error': True})
    else:
        form = UploadFileForm()
    return render(request, 'brca/home.html', {'form': form})
