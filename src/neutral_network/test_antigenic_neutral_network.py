import unittest

import numpy as np


class NeutralNetworkTests(unittest.TestCase):
    def test_make_titer_and_distance_tables(self):
        titer_table = np.array(
            [354, 0, 0, 501, 109, 0, 0],
            [640, 0, 640, 320, 80, 320, 40],
            [1163, 0, 0, 1163, 291, 0, 73],
            [1280, 1280, 0, 2560, 53, 0, 10],
            [320, 0, 0, 640, 120, 0, 80],
        )

        distance_table = np.array(
            [0, 1/0.8, 1/0.7, 1/0.6],
            [1 / 0.5, 0, 1 / 0.3, 1 / 0.2],
        )

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
