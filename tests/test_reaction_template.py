"""
Test reaction template module

Tests for bees/reaction_template.py covering:
- EC number parsing
- EC class determination
- Template determination from EC numbers
- Product inference
- Complete reaction creation

To run: pytest -v tests/test_reaction_template.py
"""

import pytest
from bees.reaction_template import (
    ECClass,
    ReactionTemplate,
    parse_ec_number,
    get_ec_class,
    determine_template_from_ec,
    infer_products,
    create_reaction_from_database
)


class TestECNumberParsing:
    """Test EC number parsing functions."""
    
    def test_parse_ec_number_with_prefix(self):
        """Test parsing EC number with 'EC ' prefix."""
        result = parse_ec_number("EC 1.2.3.4")
        assert result == (1, 2, 3, 4)
    
    def test_parse_ec_number_without_prefix(self):
        """Test parsing EC number without prefix."""
        result = parse_ec_number("2.7.1.1")
        assert result == (2, 7, 1, 1)
    
    def test_parse_ec_number_with_whitespace(self):
        """Test parsing handles whitespace."""
        result = parse_ec_number("  EC 3.1.3.1  ")
        assert result == (3, 1, 3, 1)
    
    def test_parse_ec_number_case_insensitive(self):
        """Test parsing is case-insensitive."""
        result1 = parse_ec_number("ec 1.1.1.1")
        result2 = parse_ec_number("EC 1.1.1.1")
        result3 = parse_ec_number("Ec 1.1.1.1")
        assert result1 == result2 == result3 == (1, 1, 1, 1)
    
    def test_parse_ec_number_invalid_format(self):
        """Test parsing raises error for invalid format."""
        with pytest.raises(ValueError, match="Invalid EC number format"):
            parse_ec_number("invalid")
        
        with pytest.raises(ValueError):
            parse_ec_number("1.2.3")  # Missing component
        
        with pytest.raises(ValueError):
            parse_ec_number("1.2.3.4.5")  # Too many components


class TestECClass:
    """Test EC class determination."""
    
    def test_get_ec_class_oxidoreductase(self):
        """Test EC class 1 (oxidoreductases)."""
        ec_class = get_ec_class("EC 1.1.1.1")
        assert ec_class == ECClass.OXIDOREDUCTASE
    
    def test_get_ec_class_transferase(self):
        """Test EC class 2 (transferases)."""
        ec_class = get_ec_class("EC 2.7.1.1")
        assert ec_class == ECClass.TRANSFERASE
    
    def test_get_ec_class_hydrolase(self):
        """Test EC class 3 (hydrolases)."""
        ec_class = get_ec_class("EC 3.1.3.1")
        assert ec_class == ECClass.HYDROLASE
    
    def test_get_ec_class_lyase(self):
        """Test EC class 4 (lyases)."""
        ec_class = get_ec_class("EC 4.2.1.1")
        assert ec_class == ECClass.LYASE
    
    def test_get_ec_class_isomerase(self):
        """Test EC class 5 (isomerases)."""
        ec_class = get_ec_class("EC 5.3.1.1")
        assert ec_class == ECClass.ISOMERASE
    
    def test_get_ec_class_ligase(self):
        """Test EC class 6 (ligases)."""
        ec_class = get_ec_class("EC 6.1.1.1")
        assert ec_class == ECClass.LIGASE
    
    def test_get_ec_class_translocase(self):
        """Test EC class 7 (translocases)."""
        ec_class = get_ec_class("EC 7.1.1.1")
        assert ec_class == ECClass.TRANSLOCASE
    
    def test_get_ec_class_invalid(self):
        """Test invalid EC number returns UNKNOWN."""
        ec_class = get_ec_class("invalid")
        assert ec_class == ECClass.UNKNOWN


class TestReactionTemplate:
    """Test ReactionTemplate dataclass."""
    
    def test_reaction_template_defaults(self):
        """Test default values for optional fields."""
        template = ReactionTemplate(
            template_type="generic",
            ec_class=ECClass.UNKNOWN
        )
        assert template.reactants == []
        assert template.products == []
        assert template.stoichiometry == {}
        assert template.cofactors == []
        assert template.description == ""
        assert template.reversible is True

    def test_reaction_template_repr(self):
        """Test string representation."""
        template = ReactionTemplate(
            template_type="hydrolysis",
            ec_class=ECClass.HYDROLASE,
            reactants=["ATP"],
            products=["ADP", "Pi"]
        )
        repr_str = repr(template)
        assert "hydrolysis" in repr_str
        assert "HYDROLASE" in repr_str


class TestDetermineTemplate:
    """Test template determination from EC numbers."""
    
    def test_determine_template_oxidoreductase(self):
        """Test oxidoreductase template."""
        template = determine_template_from_ec("EC 1.1.1.1")
        assert template.ec_class == ECClass.OXIDOREDUCTASE
        assert "oxidation" in template.template_type.lower()
        assert template.reversible is True
    
    def test_determine_template_transferase_kinase(self):
        """Test transferase template (kinase)."""
        template = determine_template_from_ec("EC 2.7.1.1")
        assert template.ec_class == ECClass.TRANSFERASE
        assert template.template_type == "phosphorylation"
        assert "ATP" in template.cofactors
        assert "Mg2+" in template.cofactors
        assert template.reversible is True

    def test_determine_template_hydrolase(self):
        """Test hydrolase template."""
        template = determine_template_from_ec("EC 3.1.3.1")
        assert template.ec_class == ECClass.HYDROLASE
        assert "hydrolysis" in template.template_type.lower()
        assert template.reversible is True
    
    def test_determine_template_lyase(self):
        """Test lyase template."""
        template = determine_template_from_ec("EC 4.2.1.1")
        assert template.ec_class == ECClass.LYASE
        assert template.template_type == "elimination"
        assert template.reversible is True
    
    def test_determine_template_isomerase(self):
        """Test isomerase template."""
        template = determine_template_from_ec("EC 5.3.1.1")
        assert template.ec_class == ECClass.ISOMERASE
        assert template.template_type == "isomerization"
        assert template.reversible is True
    
    def test_determine_template_ligase(self):
        """Test ligase template."""
        template = determine_template_from_ec("EC 6.1.1.1")
        assert template.ec_class == ECClass.LIGASE
        assert template.template_type == "ligation"
        assert "ATP" in template.cofactors
    
    def test_determine_template_invalid_ec(self):
        """Test template for invalid EC number."""
        template = determine_template_from_ec("invalid")
        assert template.ec_class == ECClass.UNKNOWN
        assert template.template_type == "generic"


class TestInferProducts:
    """Test product inference."""
    
    def test_infer_products_phosphorylation(self):
        """Test phosphorylation product inference."""
        template = ReactionTemplate(
            template_type="phosphorylation",
            ec_class=ECClass.TRANSFERASE
        )
        products = infer_products(
            substrate="Glucose",
            enzyme_label="Hexokinase",
            template=template,
        )
        assert any("phosphate" in p.lower() for p in products)
        assert "ADP" in products
    
    def test_infer_products_hydrolysis(self):
        """Test hydrolysis product inference."""
        template = ReactionTemplate(
            template_type="hydrolysis",
            ec_class=ECClass.HYDROLASE
        )
        products = infer_products(
            substrate="ATP",
            enzyme_label="ATPase",
            template=template
        )
        assert len(products) > 0
        assert any("hydrolyzed" in p.lower() or "ADP" in p or "Pi" in p for p in products)
    
    def test_infer_products_oxidation(self):
        """Test oxidation product inference."""
        template = ReactionTemplate(
            template_type="oxidation_CH-OH",
            ec_class=ECClass.OXIDOREDUCTASE
        )
        products = infer_products(
            substrate="Lactate",
            enzyme_label="LDH",
            template=template,
        )
        assert products == ["Lactate (oxidized)"]
    
    def test_infer_products_isomerization(self):
        """Test isomerization product inference."""
        template = ReactionTemplate(
            template_type="isomerization",
            ec_class=ECClass.ISOMERASE
        )
        products = infer_products(
            substrate="Glucose-6-phosphate",
            enzyme_label="PGI",
            template=template
        )
        assert len(products) > 0
        assert any("isomer" in p.lower() for p in products)
    
    def test_infer_products_generic(self):
        """Test generic product inference."""
        template = ReactionTemplate(
            template_type="generic",
            ec_class=ECClass.UNKNOWN
        )
        products = infer_products(
            substrate="Substrate",
            enzyme_label="Enzyme",
            template=template
        )
        assert len(products) > 0
        assert any("product" in p.lower() for p in products)


class TestCreateReactionFromDatabase:
    """Template reaction built from the EC class (the database row overwrites it when present)."""

    def test_create_reaction_inferred_products(self):
        reaction = create_reaction_from_database(
            substrate="UnknownSubstrate",
            enzyme_label="UnknownEnzyme",
            ec_number="EC 3.1.3.1",
        )
        assert reaction.reactants == ["UnknownSubstrate"]
        assert len(reaction.products) > 0
        assert reaction.stoichiometry["UnknownSubstrate"] == -1

    def test_ligase_adds_atp_reactant_and_its_products(self):
        """Ligation: ATP joins the reactants and its ADP + Pi appear once each."""
        reaction = create_reaction_from_database(
            substrate="X", enzyme_label="Ligase", ec_number="EC 6.1.1.1",
        )
        assert reaction.reactants == ["X", "ATP"]
        assert reaction.products == ["X (ligated)", "ADP", "Pi"]


class TestIntegration:
    """Integration tests."""
    
    def test_full_workflow_hexokinase(self):
        """Test complete workflow for hexokinase reaction."""
        # Parse EC number
        ec_class = get_ec_class("EC 2.7.1.1")
        assert ec_class == ECClass.TRANSFERASE
        
        # Determine template
        template = determine_template_from_ec("EC 2.7.1.1")
        assert template.template_type == "phosphorylation"
        
        # Kinase: ATP is added as a reactant; phosphorylation yields ADP but no Pi.
        reaction = create_reaction_from_database(
            substrate="Glucose",
            enzyme_label="Hexokinase",
            ec_number="EC 2.7.1.1",
        )
        assert reaction.reactants == ["Glucose", "ATP"]
        assert reaction.products == ["Glucose-phosphate", "ADP"]
        assert sum(reaction.stoichiometry.values()) == 0  # Balanced




