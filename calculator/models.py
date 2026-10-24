from django.db import models

class MoleculeQuery(models.Model):
    smiles = models.CharField(max_length=1000)
    formula = models.CharField(max_length=100, blank=True)
    molecular_weight = models.FloatField(null=True, blank=True)
    logp = models.FloatField(null=True, blank=True)
    tpsa = models.FloatField(null=True, blank=True)
    h_bond_donors = models.IntegerField(null=True, blank=True)
    h_bond_acceptors = models.IntegerField(null=True, blank=True)
    rotatable_bonds = models.IntegerField(null=True, blank=True)
    qed_score = models.FloatField(null=True, blank=True)
    lipinski_pass = models.BooleanField(default=True)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ['-created_at']

    def __str__(self):
        return f"{self.formula or 'Molecule'} - {self.smiles[:30]}"
