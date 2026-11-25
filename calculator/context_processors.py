import os
from django.conf import settings


def oauth_status(request):
    """
    Context processor that exposes OAuth configuration state to templates.
    Allows UI to gracefully handle missing/placeholder credentials without
    redirecting users to external provider 401 error pages.
    """
    google_id = getattr(settings, 'GOOGLE_CLIENT_ID', '') or ''
    github_id = getattr(settings, 'GITHUB_CLIENT_ID', '') or ''

    google_id = google_id.strip()
    github_id = github_id.strip()

    google_configured = bool(google_id and google_id != 'google-client-id-placeholder')
    github_configured = bool(github_id and github_id != 'github-client-id-placeholder')

    site_url = getattr(settings, 'SITE_URL', 'https://ilovechem.onrender.com').rstrip('/')
    
    # In local development, prefer the actual request origin for testing redirects
    if request and hasattr(request, 'get_host'):
        host = request.get_host()
        if 'localhost' in host or '127.0.0.1' in host:
            scheme = 'https' if request.is_secure() else 'http'
            site_url = f"{scheme}://{host}"

    google_callback_url = f"{site_url}/accounts/google/login/callback/"
    github_callback_url = f"{site_url}/accounts/github/login/callback/"

    return {
        'GOOGLE_AUTH_CONFIGURED': google_configured,
        'GITHUB_AUTH_CONFIGURED': github_configured,
        'GOOGLE_CALLBACK_URL': google_callback_url,
        'GITHUB_CALLBACK_URL': github_callback_url,
        'SITE_URL': site_url,
    }
