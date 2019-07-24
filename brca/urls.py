from django.conf.urls import url
from django.contrib import admin
from brca.views import home, user_login, user_logout, download

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^$', home, name='home'),
    url(r'^login/$', user_login, name='login'),
    url(r'^logout/$', user_logout, name='logout'),
    url(r'^download/$', download, name='download'),
]
