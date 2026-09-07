import unittest
import tempfile
from pathlib import Path
from mapq_calibrate import fit, upper_error, observations

class CalibrationTests(unittest.TestCase):
    def test_no_errors_does_not_imply_perfect_confidence(self):
        counts, errors = [0]*61, [0]*61
        counts[60] = 1000
        caps = fit(counts, errors)
        self.assertEqual(caps[60], 24)
        self.assertTrue(all(a <= b for a, b in zip(caps, caps[1:])))
        self.assertEqual(caps[0], 0)
        self.assertEqual(caps[59], 0)

    def test_errors_lower_confidence(self):
        self.assertGreater(upper_error(10, 100), 0.1)
        self.assertEqual(upper_error(0, 0), 1)
        counts, errors = [1000]*61, [100]*61
        self.assertTrue(all(cap <= q and cap <= 9 for q, cap in enumerate(fit(counts, errors))))

    def test_secondary_is_not_counted_as_another_read(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            truth, sam = root / 'truth', root / 'sam'
            truth.write_text('r\tchr\t10\t20\t+\n')
            sam.write_text('r\t0\tchr\t11\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n'
                           'r\t256\tother\t11\t0\t10M\t*\t0\t0\tAAAAAAAAAA\t*\n')
            counts, errors, unmapped = observations(sam, truth, 0)
            self.assertEqual((sum(counts), sum(errors), unmapped), (1, 0, 0))
            sam.write_text('')
            with self.assertRaisesRegex(ValueError, 'missing'):
                observations(sam, truth, 0)

if __name__ == '__main__':
    unittest.main()
