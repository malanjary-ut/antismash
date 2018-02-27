# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.detection.clusterfinder import probabilistic

class TestForwardBackward(unittest.TestCase):
    def test_forward_backward_short(self):
        observations = ['PF07070', 'PF01494', 'PF01266', 'PF01266', 'PF00070', 'PF03486',
                        'PF13450', 'PF01946', 'PF00890', 'PF13738', 'PF02909', 'PF00440', 'PF01494']
        expected = [0.35264287563911612, 0.9833087885291899, 0.99836086476683583, 0.99974635120762612,
                    0.9998620593169788, 0.99996163308083286, 0.99894961193916476, 0.99998899113780604,
                    0.99977313270475954, 0.99675039081490246, 0.99585473875179864, 0.99309325133832571,
                    0.99462384810384274]
        results = probabilistic.get_pfam_probabilities(observations)
        for exp, res in zip(expected, results):
            self.assertAlmostEqual(res, exp, places=6)

    def test_forward_backward_long(self):
        observations = ['PF07070', 'PF01494', 'PF01266', 'PF01266', 'PF00070', 'PF03486', 'PF13450', 'PF01946', 'PF00890', 'PF13738', 'PF02909', 'PF00440', 'PF01494', 'PF01266', 'PF13450', 'PF00890', 'PF02909', 'PF00440', 'PF01738', 'PF12695', 'PF12697', 'PF02954', 'PF14532', 'PF00158', 'PF13555', 'PF00106', 'PF13561', 'PF08659', 'PF01370', 'PF04909', 'PF01883', 'PF08240', 'PF00107', 'PF13602', 'PF01488', 'PF03992', 'PF14602', 'PF00132', 'PF00126', 'PF03466', 'PF01408', 'PF12802', 'PF01047', 'PF13463', 'PF00392', 'PF01978', 'PF09339', 'PF13412', 'PF12681', 'PF12681', 'PF00903', 'PF12833', 'PF12697', 'PF00561', 'PF12695', 'PF00975', 'PF07819', 'PF04738', 'PF14028', 'PF04737', 'PF05147', 'PF00528', 'PF00528', 'PF01547', 'PF13416', 'PF02065', 'PF01055', 'PF05691', 'PF00480', 'PF12802', 'PF09339', 'PF01978', 'PF01047', 'PF13412', 'PF01408', 'PF02894', 'PF03575', 'PF00293']
        expected = [0.35264287565878027, 0.98330878858410675, 0.99836086482622133, 0.99974635130652934, 0.99986205984404075, 0.99996163750399825, 0.99894984883917826, 0.99998945204058098, 0.99982427558092934, 0.997717586283285, 0.9977047314325167, 0.99701064972410447, 0.99988601837824154, 0.99976094552118522, 0.99847552518073068, 0.99925735569098106, 0.9941011523344635, 0.98321238809321387, 0.97911139960817117, 0.9334539912023212, 0.89142875588562287, 0.85279817889134946, 0.86982292267755557, 0.889461855725455, 0.9255595986482984, 0.96503848140239634, 0.97731600489653081, 0.99215604084380193, 0.98085602333296074, 0.45929328902331196, 0.010091418946925796, 0.056192003453458393, 0.084587568844424865, 0.097976189826937801, 0.11072136967272987, 0.078482888209672863, 0.056269642937123915, 0.033413298557164094, 0.032790576375946068, 0.063687700754878424, 0.27490446041872046, 0.33788837374664288, 0.40315265824948354, 0.46969970186839177, 0.53922188388041492, 0.88834972325210593, 0.95197541811712882, 0.94465857049585655, 0.9397400328261265, 0.93719200692148119, 0.93700009205878843, 0.94896013281974367, 0.96348126512721111, 0.98064555832182898, 0.98890088972624079, 0.99955990865971855, 0.99997878239376492, 0.99999511683006481, 0.99895421358358483, 0.99996454511323096, 0.99429232686214875, 0.19332932016966067, 0.051452598190064543, 0.027056595203090651, 0.02802568314805751, 0.027562353414535199, 0.036118414531851627, 0.2213993477782488, 0.41344942705203008, 0.60771773127165207, 0.8107132407900518, 0.81030452748214221, 0.74676442242084751, 0.68575322661921312, 0.62820362948801423, 0.43490097563133179, 0.0072909451985659996, 0.010835422528344544]
        results = probabilistic.get_pfam_probabilities(observations)
        for exp, res in zip(expected, results):
            self.assertAlmostEqual(res, exp, places=6)