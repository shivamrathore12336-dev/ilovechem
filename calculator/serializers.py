from rest_framework import serializers
from .models import MoleculeQuery

class MoleculeQuerySerializer(serializers.ModelSerializer):
    class Meta:
        model = MoleculeQuery
        fields = '__all__'
