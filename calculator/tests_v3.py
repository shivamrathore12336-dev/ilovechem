import json
from django.test import TestCase, Client
from calculator.medchem_engine import analyze_medchem_profile, batch_screen_smiles, compute_ligand_efficiency
from calculator.qm_linter import lint_qm_geometry


class V3ArchitectureTests(TestCase):
    def setUp(self):
        self.client = Client()

    def test_healthz_endpoint(self):
        resp = self.client.get("/healthz/")
        self.assertEqual(resp.status_code, 200)
        self.assertEqual(resp.json().get("status"), "ok")

    def test_pains_alert_and_boiled_egg_svg(self):
        # Clean Aspirin
        aspirin = analyze_medchem_profile("CC(=O)OC1=CC=CC=C1C(=O)O")
        self.assertIsNotNone(aspirin)
        self.assertEqual(aspirin["alerts"]["pains_count"], 0)
        self.assertIn("<svg", aspirin["boiled_egg_svg"])
        self.assertEqual(aspirin["boiled_egg"]["gi_absorption"], "High")

    def test_ligand_efficiency_calculation(self):
        le_data = compute_ligand_efficiency(affinity_kcal=-8.2, heavy_atoms=22, logp=2.4)
        self.assertTrue(le_data["le_pass"])
        self.assertGreater(le_data["le"], 0.35)

    def test_batch_screen_api(self):
        resp = self.client.post(
            "/api/medchem/screen/",
            data=json.dumps({"smiles_list": "CC(=O)OC1=CC=CC=C1C(=O)O Aspirin\nINVALID_SMILES BadMol"}),
            content_type="application/json",
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.json()
        self.assertEqual(payload["total_processed"], 1)
        self.assertEqual(payload["failed_count"], 1)

    def test_qm_linter_organometallic_heavy_atom(self):
        # Cisplatin: PtCl2(NH3)2, total Z = 78 + 2*17 + 2*7 + 6*1 = 132 electrons
        atoms = [
            ("Pt", 0.0, 0.0, 0.0),
            ("Cl", 2.3, 0.0, 0.0),
            ("Cl", 0.0, 2.3, 0.0),
            ("N", -2.0, 0.0, 0.0),
            ("N", 0.0, -2.0, 0.0),
            ("H", -2.3, 0.5, 0.5), ("H", -2.3, -0.5, 0.5), ("H", -2.3, 0.0, -0.8),
            ("H", 0.5, -2.3, 0.5), ("H", -0.5, -2.3, 0.5), ("H", 0.0, -2.3, -0.8),
        ]
        res = lint_qm_geometry(atoms, charge=0, multiplicity=1, engine="orca")
        self.assertTrue(res["valid"])
        self.assertEqual(res["electron_count"], 132)

    def test_qm_linter_catches_clash_and_parity_error(self):
        # Two atoms 0.2 A apart + invalid singlet for odd-electron radical (CH3 radical = 9e-)
        res = lint_qm_geometry(
            atoms=[("C", 0.0, 0.0, 0.0), ("H", 0.2, 0.0, 0.0), ("H", 1.1, 0.0, 0.0), ("H", 0.0, 1.1, 0.0)],
            charge=0,
            multiplicity=1,
            engine="orca",
        )
        self.assertFalse(res["valid"])
        self.assertTrue(any("Parity Mismatch" in e for e in res["errors"]))
        self.assertTrue(any("clashes" in e for e in res["errors"]))
        self.assertIn("#SBATCH", res["slurm_script"])

    def test_qm_lint_api_endpoint(self):
        resp = self.client.post(
            "/api/qm/lint/",
            data=json.dumps({
                "atoms": [["C", 0, 0, 0], ["H", 1.09, 0, 0], ["H", -0.36, 1.03, 0], ["H", -0.36, -0.51, 0.89], ["H", -0.36, -0.51, -0.89]],
                "charge": 0,
                "multiplicity": 1,
                "engine": "orca"
            }),
            content_type="application/json",
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.json()
        self.assertTrue(payload["valid"])
        self.assertEqual(payload["atom_count"], 5)

    def test_medchem_single_screen_api(self):
        resp = self.client.post(
            "/api/medchem/screen/",
            data=json.dumps({"smiles": "c1ccccc1"}),
            content_type="application/json",
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.json()
        self.assertIn("logp", payload)
        self.assertIn("boiled_egg", payload)

    def test_medchem_invalid_smiles(self):
        resp = self.client.post(
            "/api/medchem/screen/",
            data=json.dumps({"smiles": "NOT_A_MOLECULE_XYZ"}),
            content_type="application/json",
        )
        self.assertEqual(resp.status_code, 400)

    def test_properties_view_renders_boiled_egg_and_alerts(self):
        resp = self.client.get("/properties/")
        self.assertEqual(resp.status_code, 200)
        self.assertIn("BOILED-Egg Absorption Model", resp.content.decode("utf-8"))
        self.assertIn("PAINS Filter", resp.content.decode("utf-8"))

    def test_api_calculate_returns_medchem(self):
        resp = self.client.post(
            "/api/calculate/",
            data=json.dumps({"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}),
            content_type="application/json",
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertIn("medchem", data)
        self.assertEqual(data["medchem"]["alerts"]["pains_count"], 0)

    def test_analyze_none_and_empty(self):
        self.assertIsNone(analyze_medchem_profile(None))
        self.assertIsNone(analyze_medchem_profile(""))
        self.assertIsNone(analyze_medchem_profile("   "))
