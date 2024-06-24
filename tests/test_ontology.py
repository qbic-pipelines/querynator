#!/usr/bin/env python

"""Tests for Ontology classes."""

import unittest

from querynator.helper_functions import ontology


class testDiseaseOntology(unittest.TestCase):
    """Test the DiseaseOntology class"""

    def setUp(self):
        """Set up test fixtures"""
        self.doid = ontology.DiseaseOntology("querynator/helper_functions/doid.obo")

    def test_contains(self):
        """Test term containment in DiseaseOntology"""
        self.assertIn("DOID:4", self.doid.ids)
        self.assertIn("DOID:4", self.doid.terms)
        self.assertIn("DOID:4", self.doid)
        self.assertIn("DOID:4947", self.doid)

    def test_get(self):
        """Test term retrieval in DiseaseOntology"""
        term = self.doid.get("DOID:4947")
        self.assertEqual(term, self.doid.get(4947))
        self.assertEqual(term, self.doid.get("cholangiocarcinoma"))
        self.assertEqual(term, self.doid.get("Cholangiocarcinoma"))

        self.assertIsNone(self.doid.get("DOID:123456"))
        self.assertIsNone(self.doid.get_from_name("foo"))

    def test_relationships(self):
        """Test relationships in DiseaseOntology"""
        self.assertIn(
            self.doid.get("bile duct adenocarcinoma"), self.doid.get("cholangiocarcinoma").get_related_terms()
        )
        self.assertIn(self.doid.get("bile duct adenocarcinoma"), self.doid.get_all_ancestors("cholangiocarcinoma"))
        self.assertIn(
            self.doid.get("cholangiocarcinoma"), self.doid.get_all_ancestors("cholangiocarcinoma", includeSelf=True)
        )

        self.assertNotIn(
            self.doid.get("cholangiocarcinoma"), self.doid.get_all_ancestors("cholangiocarcinoma", includeSelf=False)
        )
        self.assertNotIn(
            self.doid.get("cholangiocarcinoma"),
            self.doid.get_all_ancestors("bile duct adenocarcinoma", includeSelf=True),
        )


if __name__ == "__main__":
    unittest.main()
