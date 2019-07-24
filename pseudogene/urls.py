from django.conf.urls import url
from django.contrib import admin

from pseudogene.views import search_form, test_forms

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^$', search_form, name='index'),
    url(r'^test/$', test_forms, name='test'),
]
