from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path('admin/', admin.site.urls),
    path('accounts/', include('allauth.urls')),
    # Every public route (pages, robots.txt, sitemap.xml, og-image, APIs) lives in calculator/urls.py.
    path('', include('calculator.urls')),
]
