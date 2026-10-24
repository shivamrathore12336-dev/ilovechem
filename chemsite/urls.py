from django.contrib import admin
from django.urls import path, include
from calculator import studio_views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('robots.txt', studio_views.robots_txt_view, name='robots_txt'),
    path('sitemap.xml', studio_views.sitemap_xml_view, name='sitemap_xml'),
    path('', include('calculator.urls')),
]
