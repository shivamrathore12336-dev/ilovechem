from django.db import migrations


class Migration(migrations.Migration):
    """Drop the unused MoleculeQuery history table to honor the zero-storage policy."""

    dependencies = [
        ('calculator', '0001_initial'),
    ]

    operations = [
        migrations.DeleteModel(
            name='MoleculeQuery',
        ),
    ]
