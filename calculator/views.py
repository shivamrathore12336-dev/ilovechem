from django.shortcuts import render
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from .chem_utils import (
    validate_smiles,
    calculate_properties,
    generate_3d_molblock,
    generate_2d_base64_image
)
from .models import MoleculeQuery
from .serializers import MoleculeQuerySerializer

def home_view(request):
    """
    Main interactive dashboard view.
    Handles form submission, computes properties and 3D structure,
    and returns context for instant rendering.
    """
    context = {}
    default_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    selected_smiles = default_smiles

    if request.method == "POST":
        selected_smiles = request.POST.get("smiles", "").strip()

    if selected_smiles:
        context["input_smiles"] = selected_smiles
        if not validate_smiles(selected_smiles):
            context["error"] = f"Invalid SMILES string: '{selected_smiles}'. Please verify the chemical notation."
        else:
            props = calculate_properties(selected_smiles)
            if "error" in props:
                context["error"] = props["error"]
            else:
                context["properties"] = props
                context["img_2d"] = generate_2d_base64_image(selected_smiles)
                structure_res = generate_3d_molblock(selected_smiles)
                if "mol_block" in structure_res:
                    context["mol_block"] = structure_res["mol_block"]
                else:
                    context["structure_warning"] = structure_res.get("error", "3D generation unavailable")

                # Store history asynchronously or directly in SQLite
                try:
                    MoleculeQuery.objects.create(
                        smiles=selected_smiles,
                        formula=props.get("formula", ""),
                        molecular_weight=props.get("molecular_weight"),
                        logp=props.get("logp"),
                        tpsa=props.get("tpsa"),
                        h_bond_donors=props.get("h_bond_donors"),
                        h_bond_acceptors=props.get("h_bond_acceptors"),
                        rotatable_bonds=props.get("rotatable_bonds"),
                        qed_score=props.get("qed_score"),
                        lipinski_pass=props.get("lipinski_pass", True)
                    )
                except Exception:
                    pass

    context["recent_queries"] = MoleculeQuery.objects.all()[:6]
    return render(request, "calculator/index.html", context)


@api_view(['POST'])
def api_calculate_properties(request):
    """
    JSON API Endpoint for external tools and programmatic use.
    POST Body: {"smiles": "..."}
    """
    smiles = request.data.get("smiles", "").strip()
    if not smiles or not validate_smiles(smiles):
        return Response(
            {"error": "Valid 'smiles' string is required."},
            status=status.HTTP_400_BAD_REQUEST
        )

    props = calculate_properties(smiles)
    if "error" in props:
        return Response(props, status=status.HTTP_400_BAD_REQUEST)

    structure_3d = generate_3d_molblock(smiles)
    img_2d = generate_2d_base64_image(smiles)

    return Response({
        "properties": props,
        "mol_block": structure_3d.get("mol_block", ""),
        "image_2d": img_2d
    }, status=status.HTTP_200_OK)


@api_view(['GET'])
def api_history(request):
    """Returns recent queries."""
    queries = MoleculeQuery.objects.all()[:20]
    serializer = MoleculeQuerySerializer(queries, many=True)
    return Response(serializer.data)
