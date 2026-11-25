from django.test import TestCase, Client, override_settings
from django.contrib.auth.models import User
from calculator.context_processors import oauth_status
from django.test.client import RequestFactory


class AuthenticationSystemTests(TestCase):
    def setUp(self):
        self.client = Client()
        self.rf = RequestFactory()

    def test_auth_login_page_loads(self):
        resp = self.client.get('/accounts/login/')
        self.assertEqual(resp.status_code, 200)

    def test_anonymous_user_sees_sign_in_button_in_header(self):
        resp = self.client.get('/')
        self.assertEqual(resp.status_code, 200)
        content = resp.content.decode('utf-8')
        self.assertIn('Sign In', content)
        self.assertIn('Sign In to iLoveChem', content)
        self.assertIn('Continue with Google', content)
        self.assertIn('Continue with GitHub', content)
        self.assertIn('Continue as Guest Researcher', content)

    def test_authenticated_user_sees_avatar_and_sign_out(self):
        user = User.objects.create_user(username='testchemist', email='chemist@university.edu', password='password123')
        self.client.force_login(user)
        resp = self.client.get('/')
        self.assertEqual(resp.status_code, 200)
        content = resp.content.decode('utf-8')
        self.assertIn('chemist@university.edu', content)
        self.assertIn('/accounts/logout/', content)

    def test_google_and_github_login_endpoints_exist(self):
        resp_google = self.client.get('/accounts/google/login/')
        self.assertIn(resp_google.status_code, [200, 302])

        resp_github = self.client.get('/accounts/github/login/')
        self.assertIn(resp_github.status_code, [200, 302])

    def test_logout_redirects_cleanly(self):
        user = User.objects.create_user(username='logoutuser', email='logout@chem.org', password='password123')
        self.client.force_login(user)
        resp = self.client.get('/accounts/logout/')
        self.assertEqual(resp.status_code, 302)
        self.assertEqual(resp.url, '/')

    def test_guest_login_creates_user_and_session(self):
        resp = self.client.post('/accounts/guest-login/', {'next': '/docking/'})
        self.assertEqual(resp.status_code, 302)
        self.assertEqual(resp.url, '/docking/')

        # Verify session is authenticated
        resp_home = self.client.get('/')
        content = resp_home.content.decode('utf-8')
        self.assertIn('guest.researcher@ilovechem.org', content)

    def test_guest_login_safe_next_redirect(self):
        # Malicious or external next parameter should fall back to '/'
        resp = self.client.post('/accounts/guest-login/', {'next': 'https://evil.com'})
        self.assertEqual(resp.status_code, 302)
        self.assertEqual(resp.url, '/')

    def test_oauth_status_context_processor(self):
        req = self.rf.get('/')
        ctx = oauth_status(req)
        self.assertIn('GOOGLE_AUTH_CONFIGURED', ctx)
        self.assertIn('GITHUB_AUTH_CONFIGURED', ctx)
        self.assertIn('GOOGLE_CALLBACK_URL', ctx)
        self.assertTrue(ctx['GOOGLE_CALLBACK_URL'].endswith('/accounts/google/login/callback/'))

    @override_settings(GOOGLE_CLIENT_ID='actual-google-client-id.apps.googleusercontent.com')
    def test_google_auth_configured_when_valid_client_id_present(self):
        req = self.rf.get('/')
        ctx = oauth_status(req)
        self.assertTrue(ctx['GOOGLE_AUTH_CONFIGURED'])
