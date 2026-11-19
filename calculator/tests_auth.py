from django.test import TestCase, Client
from django.contrib.auth.models import User


class AuthenticationSystemTests(TestCase):
    def setUp(self):
        self.client = Client()

    def test_auth_login_page_loads(self):
        resp = self.client.get('/accounts/login/')
        self.assertEqual(resp.status_code, 200)

    def test_anonymous_user_sees_sign_in_button_in_header(self):
        resp = self.client.get('/')
        self.assertEqual(resp.status_code, 200)
        self.assertIn('Sign In', resp.content.decode('utf-8'))
        self.assertIn('Sign In to iLoveChem', resp.content.decode('utf-8'))
        self.assertIn('Continue with Google', resp.content.decode('utf-8'))
        self.assertIn('Continue with GitHub', resp.content.decode('utf-8'))

    def test_authenticated_user_sees_avatar_and_sign_out(self):
        user = User.objects.create_user(username='testchemist', email='chemist@university.edu', password='password123')
        self.client.force_login(user)
        resp = self.client.get('/')
        self.assertEqual(resp.status_code, 200)
        content = resp.content.decode('utf-8')
        self.assertIn('chemist@university.edu', content)
        self.assertIn('/accounts/logout/', content)

    def test_google_and_github_login_endpoints_exist(self):
        # GET request to google login triggers redirect or page
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
