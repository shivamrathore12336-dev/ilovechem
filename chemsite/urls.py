from django.contrib import admin
from django.urls import path, include
from calculator import studio_views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('robots.txt', studio_views.robots_txt_view, name='robots_txt'),
    path('sitemap.xml', studio_views.sitemap_xml_view, name='sitemap_xml'),
    path('google889b369925c8069a.html', studio_views.google_verification_view, name='google_verify'),
    path('og-image.png', studio_views.og_image_view, name='og_image'),
    path('downloads/', studio_views.downloads_view, name='downloads_hub'),
    path('', include('calculator.urls')),
]
