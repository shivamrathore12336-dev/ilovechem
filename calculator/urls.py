from django.urls import path
from . import views
from . import studio_views

urlpatterns = [
    # Universal Converter Studio (Main Hub)
    path('', studio_views.studio_view, name='studio_home'),
    path('downloads/', studio_views.downloads_view, name='downloads_page'),
    path('robots.txt', studio_views.robots_txt_view, name='robots_txt'),
    path('sitemap.xml', studio_views.sitemap_xml_view, name='sitemap_xml'),
    path('google889b369925c8069a.html', studio_views.google_verification_view, name='google_verify_calc'),
    path('api/convert/', studio_views.api_convert_file, name='api_convert'),
    path('audit/bundle/', studio_views.audit_bundle_view, name='audit_bundle'),
    
    # Original Descriptors & API
    path('properties/', views.home_view, name='properties_home'),
    path('api/calculate/', views.api_calculate_properties, name='api_calculate'),
    path('api/history/', views.api_history, name='api_history'),
]
