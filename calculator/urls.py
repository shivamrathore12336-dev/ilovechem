from django.urls import path
from . import views
from . import studio_views
from . import docking_views

urlpatterns = [
    # Universal Converter Studio (Main Hub)
    path('', studio_views.studio_view, name='studio_home'),
    path('downloads/', studio_views.downloads_view, name='downloads_page'),
    path('about/', studio_views.about_view, name='about_page'),
    path('robots.txt', studio_views.robots_txt_view, name='robots_txt'),
    path('sitemap.xml', studio_views.sitemap_xml_view, name='sitemap_xml'),
    path('google889b369925c8069a.html', studio_views.google_verification_view, name='google_verify_calc'),
    path('api/convert/', studio_views.api_convert_file, name='api_convert'),
    # Original Descriptors & API
    path('properties/', views.home_view, name='properties_home'),
    path('api/calculate/', views.api_calculate_properties, name='api_calculate'),
    # Inbuilt Molecular Docking Studio (Phase 2)
    path('docking/', docking_views.docking_studio_view, name='docking_studio'),
    path('api/dock/fetch-pdb/<str:pdb_id>/', docking_views.api_fetch_rcsb_pdb, name='api_fetch_pdb'),
    path('api/dock/rcsb/<str:pdb_id>/', docking_views.api_fetch_rcsb_pdb, name='api_fetch_rcsb_pdb_alias'),
    path('api/dock/run/', docking_views.api_run_docking, name='api_run_docking'),
]
