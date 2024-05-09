import re
from abc import ABC, abstractmethod
from typing import List


class Ontology(ABC):
    """Abstract class for ontologies.

    Ontologies represent hierarchical relationships between terms.
    """

    @abstractmethod
    def __init__(self, path):
        self.terms = {}
        self.relationships = {}
        self.path = path
        self.ids = set()
        self.names = set()

    @abstractmethod
    def get_all_ancestors(self, query_term):
        pass

    @abstractmethod
    def get(self, term_id):
        pass

    def __contains__(self, item):
        return item in self.terms

    def __repr__(self):
        return f"<Ontology.{self.__class__.__name__} {self.path}>"

    class Term:
        """represents an ontology term"""

        def __init__(self, fields: dict, relationships: List["Ontology.Relationship"] = None):
            try:
                self.id = fields["id"]
            except KeyError:
                raise ValueError("The 'id' field is required")

            self.fields = fields
            self.relationships = relationships

        def __eq__(self, other):
            if isinstance(other, self.__class__):
                return self.id == other.id
            elif isinstance(other, str):
                return self.id.lower() == other.lower() or self.get("name").lower() == other.lower()
            elif isinstance(other, int):
                return self.id.split(":")[1] == other
            else:
                return False

        def __repr__(self):
            return f"<{self.id} {self.get('name')}>"

        def __hash__(self):
            return hash(self.id)

        def __getattr__(self, attr):
            return self.fields.get(attr)

        def get(self, field):
            return self.fields.get(field)

        def get_related_terms(self):
            """return the terms directly related to this term"""
            return [rel.T2 for rel in self.relationships]

    class Relationship:
        """represents a relationship between two ontology terms"""

        def __init__(self, term1, rtype, term2, id=None):
            self.id = id
            self.T1 = term1
            self.rtype = rtype
            self.T2 = term2

        def __repr__(self):
            return f"<Ontology.Relationship: {self.T1} {self.rtype} {self.T2}>"

        def __eq__(self, other):
            if isinstance(other, Ontology.Term):
                return self.T2 == other
            else:
                return hash(self) == hash(other)

        def __contains__(self, item):
            return item == self.T1 or item == self.T2

        def __hash__(self):
            return hash(f"{self.T1} {self.rtype} {self.T2}")


class DiseaseOntology(Ontology):
    """This class represents the Disease Ontology (DO)

    www.disease-ontology.org
    """

    def __init__(self, path):
        self.terms = {}
        self.ids = set()
        self.names = set()
        self.path = path

        with open(path, "r") as file:
            content = file.read()
            terms = re.findall(r"\[Term\]\n(.*?)\n\n", content, re.DOTALL)
            for term in terms:
                term_dict = {}
                lines = term.split("\n")
                for line in lines:
                    if line:
                        # parse fields into dict
                        key, value = line.split(": ", 1)
                        if key in term_dict:
                            if isinstance(term_dict[key], list):
                                term_dict[key].append(value)
                            else:
                                term_dict[key] = [term_dict[key], value]
                        else:
                            term_dict[key] = value
                # create term object
                is_a_values = term_dict.get("is_a")
                if isinstance(is_a_values, list):
                    relationships = [
                        self.Relationship(term_dict.get("id"), "is_a", other.split(" ! ")[0]) for other in is_a_values
                    ]
                elif isinstance(is_a_values, str):
                    relationships = [self.Relationship(term_dict.get("id"), "is_a", is_a_values.split(" ! ")[0])]
                else:
                    relationships = []
                new_term = self.Term(term_dict, relationships)
                # add term to the ontology
                self.terms[new_term.id] = new_term
                self.ids.add(new_term.id)
                self.names.add(new_term.get("name"))

    def get(self, query):
        """retrieve a term from the ontology

        query can be one of id (str), id (int), or name (str)
        examples: "DOID:4947", 4947, "cholangiocarcinoma"
        """
        if isinstance(query, int):
            return self.terms.get(f"DOID:{query}")
        elif isinstance(query, str):
            if query.startswith("DOID:"):
                return self.terms.get(query)
            else:
                try:
                    return self.get(int(query))
                except ValueError:
                    return self.get_from_name(query.lower())
        else:
            raise ValueError(f"Invalid query type: {type(query)}")

    def get_from_name(self, name) -> Ontology.Term:
        """retrieve a term from the ontology by name

        :param name: the name of the term to retrieve
        :type  name: str
        :return: the term object if found, None otherwise
        """
        for id, term in self.terms.items():
            if term == name:
                return term
        else:
            return None

    def get_all_ancestors(self, query_term, includeSelf=False) -> list:
        # Get the term corresponding to the query DOID
        term = self.get(query_term)
        ancestors = set([term]) if includeSelf else set()
        if term:
            visited = set()
            # Traverse the graph to find all ancestors
            stack = [term]
            while stack:
                current_term = stack.pop()
                visited.add(current_term.id)

                for parent in [self.get(rel.T2) for rel in term.relationships if rel.rtype == "is_a"]:
                    if parent.id not in visited:
                        ancestors.add(parent)
                        # Add the parent term to the stack for further traversal
                        stack.append(parent)
            return list(ancestors)
        else:
            return []

    def get_all_ancestor_names(self, query_term) -> list:
        return [self.get(ancestor).name for ancestor in self.get_all_ancestors(query_term)]
