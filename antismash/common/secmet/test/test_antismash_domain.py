# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet.feature import AntismashDomain, FeatureLocation


class TestConversion(unittest.TestCase):
    def test_conversion(self):
        domain = AntismashDomain(FeatureLocation(1, 3, 1))
        domain.domain_subtype = "subtest"
        domain.specificity = ["a", "c", "f"]
        domain.asf.add("first")
        domain.asf.add("second")

        bio = domain.to_biopython()
        assert len(bio) == 1
        new_domain = AntismashDomain.from_biopython(bio[0])
        assert new_domain.domain_subtype == domain.domain_subtype == "subtest"
        assert new_domain.specificity == domain.specificity == ["a", "c", "f"]
        assert new_domain.asf.hits == domain.asf.hits
        assert new_domain.asf.hits == ["first", "second"]