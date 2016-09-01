from django.conf.urls import url
from django.contrib import admin

from pseudogene.views import search_form

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^$', search_form, name='index'),
]
