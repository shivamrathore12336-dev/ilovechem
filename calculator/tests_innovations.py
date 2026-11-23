import json
from django.test import TestCase, Client
from calculator.lead_optimizer import generate_bioisostere_analogs
from calculator.similarity_engine import compare_library_to_reference
from calculator.pubchem_resolver import resolve_compound_identifier


class InnovativeFeatureTests(TestCase):
    def setUp(self):
        self.client = Client()

    def test_bioisostere_optimizer_proposes_tetrazole_and_amide_for_aspirin(self):
        res = generate_bioisostere_analogs("CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertIsNotNone(res)
        self.assertGreater(res["analog_count"], 0)
        rules = {a["rule"] for a in res["analogs"]}
        self.assertIn("Carboxylic Acid -> Tetrazole", rules)
        # Check delta values exist
        self.assertIn("delta_qed", res["analogs"][0])
        self.assertIn("delta_logp", res["analogs"][0])

    def test_tanimoto_similarity_ranks_identical_molecule_first(self):
        ref = "CC(=O)OC1=CC=CC=C1C(=O)O"
        res = compare_library_to_reference(ref, [
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C Caffeine",
            "CC(=O)OC1=CC=CC=C1C(=O)O Aspirin",
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O Ibuprofen",
        ])
        self.assertIsNotNone(res)
        self.assertEqual(res["results"][0]["name"], "Aspirin")
        self.assertAlmostEqual(res["results"][0]["tanimoto_ecfp4"], 1.0, places=3)
        self.assertTrue(res["results"][0]["same_murcko_scaffold"])

    def test_api_lead_optimize_endpoint(self):
        resp = self.client.post(
            "/api/medchem/optimize/",
            data=json.dumps({"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}),
            content_type="application/json"
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertIn("analogs", data)
        self.assertGreater(data["analog_count"], 0)

    def test_api_similarity_compare_endpoint(self):
        resp = self.client.post(
            "/api/medchem/similarity/",
            data=json.dumps({
                "reference_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "candidates": "CC(=O)OC1=CC=CC=C1C(=O)O Aspirin\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C Caffeine"
            }),
            content_type="application/json"
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertEqual(data["compared_count"], 2)
        self.assertEqual(data["results"][0]["name"], "Aspirin")

    def test_pubchem_resolver_handles_invalid_and_empty(self):
        self.assertIsNone(resolve_compound_identifier(None))
        self.assertIsNone(resolve_compound_identifier(""))
        self.assertIsNone(resolve_compound_identifier("   "))

    def test_bioavailability_radar_svg_generation_and_esol_logs(self):
        from calculator.medchem_engine import analyze_medchem_profile, compute_esol_logs, generate_bioavailability_radar_svg
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
        logs = compute_esol_logs(mol, 180.16, 1.31, 3)
        self.assertIsInstance(logs, float)
        self.assertLess(logs, 0.0)

        svg = generate_bioavailability_radar_svg(180.16, 1.31, 63.6, 3, 0.11, logs, label="Aspirin")
        self.assertIn("<svg", svg)
        self.assertIn("Bioavailability Radar", svg)
        self.assertIn("LIPO", svg)
        self.assertIn("POLAR", svg)

        prof = analyze_medchem_profile("CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertIn("bioavailability_radar_svg", prof)
        self.assertIn("esol_logs", prof)

    def test_home_view_renders_radar_and_bioisosteres(self):
        resp = self.client.get("/properties/?smiles=CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertEqual(resp.status_code, 200)
        content = resp.content.decode("utf-8")
        self.assertIn("Bioavailability Radar", content)
        self.assertIn("Lead Optimization &amp; Bioisosteric Replacements", content)
        self.assertIn("Carboxylic Acid -&gt; Tetrazole", content)

    def test_sitemap_and_robots_cache_control_and_lastmod(self):
        resp_sitemap = self.client.get("/sitemap.xml")
        self.assertEqual(resp_sitemap.status_code, 200)
        self.assertEqual(resp_sitemap["Cache-Control"], "public, max-age=86400")
        self.assertIn("<lastmod>", resp_sitemap.content.decode("utf-8"))

        resp_robots = self.client.get("/robots.txt")
        self.assertEqual(resp_robots.status_code, 200)
        self.assertEqual(resp_robots["Cache-Control"], "public, max-age=86400")

    def test_pages_have_canonical_and_json_ld_schemas(self):
        resp_prop = self.client.get("/properties/")
        self.assertEqual(resp_prop.status_code, 200)
        c_prop = resp_prop.content.decode("utf-8")
        self.assertIn('<link rel="canonical" href="https://ilovechem.onrender.com/properties/">', c_prop)
        self.assertIn('"@type": "FAQPage"', c_prop)

        resp_dock = self.client.get("/docking/")
        self.assertEqual(resp_dock.status_code, 200)
        c_dock = resp_dock.content.decode("utf-8")
        self.assertIn('<link rel="canonical" href="https://ilovechem.onrender.com/docking/">', c_dock)
        self.assertIn('"@type": "HowTo"', c_dock)
